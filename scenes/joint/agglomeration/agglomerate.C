// Select aligned Cartesian sibling octets and merge with OpenFOAM removeFaces.
#include "fvCFD.H"
#include "argList.H"
#include "Time.H"
#include "labelIOList.H"
#include "polyTopoChange.H"
#include "removeFaces.H"

#include <array>
#include <map>

using namespace Foam;

namespace
{
const point origin(-600, -1200, 0);
constexpr scalar tolerance = 1e-7;

bool overlaps(const point& low, const point& high, const point& a, const point& b)
{
    return low.x() < b.x() && high.x() > a.x()
        && low.y() < b.y() && high.y() > a.y()
        && low.z() < b.z() && high.z() > a.z();
}

bool cartesianCell
(
    const labelList& vertices,
    const pointField& points,
    const boolList& patchPoint,
    const point& low,
    const scalar width
)
{
    if (vertices.size() != 8) return false;
    unsigned corners = 0;
    for (const label vertex : vertices)
    {
        if (patchPoint[vertex]) return false;
        unsigned corner = 0;
        for (direction axis = 0; axis < 3; ++axis)
        {
            const scalar t = (points[vertex][axis] - low[axis]) / width;
            if (mag(t - 1) < tolerance) corner |= 1u << axis;
            else if (mag(t) >= tolerance) return false;
        }
        corners |= 1u << corner;
    }
    return corners == 255;
}

void verifyPatches
(
    const fvMesh& mesh,
    const polyBoundaryMesh& patches,
    const wordList& names,
    const wordList& types,
    const labelList& sizes,
    const faceList& oldFaces,
    const pointField& oldPoints,
    const label boundaryStart,
    const mapPolyMesh& mapping
)
{
    if (patches.size() != names.size())
        FatalErrorInFunction << "Boundary patch count changed" << exit(FatalError);
    label offset = 0;
    forAll(patches, patchi)
    {
        const polyPatch& patch = patches[patchi];
        if (patch.name() != names[patchi] || patch.type() != types[patchi]
            || patch.size() != sizes[patchi])
            FatalErrorInFunction << "Boundary patch changed" << exit(FatalError);
        forAll(patch, i)
        {
            const label facei = patch.start() + i;
            const face& before = oldFaces[offset + i];
            const face& after = mesh.faces()[facei];
            if (mapping.faceMap()[facei] != boundaryStart + offset + i
                || before.size() != after.size())
                FatalErrorInFunction << "Boundary face changed" << exit(FatalError);
            forAll(after, vertex)
            {
                const point& a = mesh.points()[after[vertex]];
                const point& b = oldPoints[before[vertex]];
                if (a.x() != b.x() || a.y() != b.y() || a.z() != b.z())
                    FatalErrorInFunction << "Boundary point changed" << exit(FatalError);
            }
        }
        offset += patch.size();
    }
}
}

int main(int argc, char* argv[])
{
    argList::noFunctionObjects();
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    const word instance = mesh.pointsInstance();
    labelIOList levels
    (
        IOobject("cellLevel", instance, polyMesh::meshSubDir, mesh, IOobject::MUST_READ)
    );
    labelIOList pointLevels
    (
        IOobject("pointLevel", instance, polyMesh::meshSubDir, mesh, IOobject::MUST_READ)
    );
    const label oldCellCount = mesh.nCells();
    if (levels.size() != oldCellCount || pointLevels.size() != mesh.nPoints())
        FatalErrorInFunction << "Invalid input refinement level size" << exit(FatalError);

    const pointField oldPoints(mesh.points());
    const polyBoundaryMesh& patches = mesh.boundaryMesh();
    const label boundaryStart = mesh.nInternalFaces();
    faceList oldFaces(mesh.nFaces() - boundaryStart);
    boolList patchPoint(mesh.nPoints(), false);
    forAll(oldFaces, i)
    {
        oldFaces[i] = mesh.faces()[boundaryStart + i];
        for (const label vertex : oldFaces[i]) patchPoint[vertex] = true;
    }
    wordList patchNames(patches.size()), patchTypes(patches.size());
    labelList patchSizes(patches.size());
    forAll(patches, i)
    {
        patchNames[i] = patches[i].name();
        patchTypes[i] = patches[i].type();
        patchSizes[i] = patches[i].size();
    }

    using Key = std::array<label, 4>;
    std::map<Key, DynamicList<label>> groups;
    const labelListList& vertices = mesh.cellPoints();
    const labelListList& neighbours = mesh.cellCells();
    const scalarField oldVolumes(mesh.V());
    const vectorField& centres = mesh.cellCentres();
    forAll(levels, celli)
    {
        const label level = levels[celli];
        if (level < 0 || level > 3) continue;
        const scalar width = std::ldexp(20., -level);
        if (mag(oldVolumes[celli] - pow3(width)) > tolerance * pow3(width)) continue;
        const vector grid = (centres[celli] - origin) / width;
        const vector index(floor(grid.x()), floor(grid.y()), floor(grid.z()));
        const point low = origin + width * vector(index.x(), index.y(), index.z());
        if (!cartesianCell(vertices[celli], oldPoints, patchPoint, low, width)) continue;
        const Key key{level, label(index.x())/2, label(index.y())/2, label(index.z())/2};
        const point parentLow = origin + 2*width * vector(key[1], key[2], key[3]);
        const point parentHigh = parentLow + 2*width * vector(1,1,1);
        // Preserve the building district, rotor and sampled wake/flight corridor.
        if (overlaps(parentLow, parentHigh, point(-270,-270,0), point(270,270,140))
            || overlaps(parentLow, parentHigh, point(-65,-480,0), point(65,-195,105)))
            continue;
        bool touchesFiner = false;
        for (const label adjacent : neighbours[celli])
            if (levels[adjacent] > level) touchesFiner = true;
        if (!touchesFiner) groups[key].append(celli);
    }

    labelList groupForCell(oldCellCount, -1);
    DynamicList<point> expectedCentres;
    DynamicList<scalar> expectedVolumes;
    labelList countByLevel(4, 0);
    for (const auto& entry : groups)
    {
        const Key& key = entry.first;
        const DynamicList<label>& cells = entry.second;
        if (cells.size() != 8) continue;
        const scalar width = std::ldexp(20., -key[0]);
        const label group = expectedCentres.size();
        expectedCentres.append(origin + width * vector(2*key[1]+1, 2*key[2]+1, 2*key[3]+1));
        expectedVolumes.append(pow3(2*width));
        ++countByLevel[key[0]];
        for (const label celli : cells) groupForCell[celli] = group;
    }
    const label nGroups = expectedCentres.size();
    Info << "Sibling octets by original level (0 through 3): " << countByLevel << nl;
    if (!nGroups)
        FatalErrorInFunction << "No eligible cell octets" << exit(FatalError);

    DynamicList<label> candidates;
    labelList countByGroup(nGroups, 0);
    for (label facei = 0; facei < boundaryStart; ++facei)
    {
        const label group = groupForCell[mesh.faceOwner()[facei]];
        if (group >= 0 && group == groupForCell[mesh.faceNeighbour()[facei]])
        {
            candidates.append(facei);
            ++countByGroup[group];
        }
    }
    for (const label nFaces : countByGroup)
        if (nFaces != 12)
            FatalErrorInFunction << "Octet does not share twelve internal faces" << exit(FatalError);

    removeFaces remover(mesh, 0.999999);
    labelList regions, masters, compatible;
    const label nRegions = remover.compatibleRemoves(candidates, regions, masters, compatible);
    if (compatible != candidates || nRegions != nGroups)
        FatalErrorInFunction << "Unexpected compatible face merges" << exit(FatalError);
    polyTopoChange topology(mesh);
    remover.setRefinement(compatible, regions, masters, topology);
    autoPtr<mapPolyMesh> mapping = topology.changeMesh(mesh, false);
    mesh.updateMesh(mapping());
    if (mapping().hasMotionPoints()) mesh.movePoints(mapping().preMotionPoints());

    verifyPatches(mesh, mesh.boundaryMesh(), patchNames, patchTypes, patchSizes,
        oldFaces, oldPoints, boundaryStart, mapping());
    if (mesh.nCells() != oldCellCount - 7*nGroups)
        FatalErrorInFunction << "Unexpected resulting cell count" << exit(FatalError);

    const scalarField& volumes = mesh.V();
    labelList remapped(mesh.nCells());
    boolList visits(nGroups, false);
    forAll(remapped, celli)
    {
        const label oldCell = mapping().cellMap()[celli];
        if (oldCell < 0 || oldCell >= oldCellCount)
            FatalErrorInFunction << "Result cell has no source cell" << exit(FatalError);
        const label group = groupForCell[oldCell];
        remapped[celli] = levels[oldCell] - (group >= 0);
        const scalar targetVolume = group >= 0 ? expectedVolumes[group] : oldVolumes[oldCell];
        if (mag(volumes[celli] - targetVolume) > tolerance * max(1., mag(targetVolume)))
            FatalErrorInFunction << "Map disagrees with geometry at cell " << celli
                << ": expected volume " << targetVolume << ", got " << volumes[celli]
                << exit(FatalError);
        if (group >= 0)
        {
            if (visits[group] || mag(mesh.cellCentres()[celli] - expectedCentres[group]) > tolerance)
                FatalErrorInFunction << "Mapped parent centre mismatch" << exit(FatalError);
            visits[group] = true;
        }
    }
    for (const bool visited : visits)
        if (!visited) FatalErrorInFunction << "Missing parent cell" << exit(FatalError);

    labelList remappedPoints(mesh.nPoints());
    forAll(remappedPoints, pointi)
    {
        const label oldPoint = mapping().pointMap()[pointi];
        if (oldPoint < 0 || oldPoint >= pointLevels.size())
            FatalErrorInFunction << "Result point has no source point" << exit(FatalError);
        remappedPoints[pointi] = pointLevels[oldPoint];
    }
    mesh.setInstance(instance);
    levels.transfer(remapped);
    pointLevels.transfer(remappedPoints);
    if (!mesh.write() || !levels.write() || !pointLevels.write())
        FatalErrorInFunction << "Failed to write mesh and mapped levels" << exit(FatalError);
    Info << "Verified unchanged patch faces and remapped cell/point levels; cells "
        << oldCellCount << " -> " << mesh.nCells() << nl;
}
