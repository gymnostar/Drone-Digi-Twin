// Agglomerate complete snapped sibling octets without changing surviving face polygons.
#include "fvCFD.H"
#include "argList.H"
#include "Time.H"
#include "labelIOList.H"
#include "polyTopoChange.H"

#include <array>
#include <map>

using namespace Foam;

namespace
{
const point origin(-600, -1200, 0);
constexpr scalar tolerance = 1e-7;
constexpr scalar allowedDisplacement = 0.18;

struct Octet
{
    std::array<label, 8> children;
    bool valid = true;
    Octet() { children.fill(-1); }
};

struct Moments
{
    scalar volume = 0;
    vector first = vector::zero;
};

struct Parent
{
    scalar volume;
    point centroid;
    Moments exact;
};

// Integrate identical oriented face fans relative to a common local origin.
void addCellMoments(Moments& total, const point& origin, const cell& faces,
    const faceList& polygons, const pointField& points, const labelList& owners,
    const label celli)
{
    for (const label facei : faces)
    {
        const face& polygon = polygons[facei];
        vector centre = vector::zero;
        for (const label vertex : polygon) centre += points[vertex] - origin;
        centre /= polygon.size();
        const scalar orientation = owners[facei] == celli ? 1 : -1;
        forAll(polygon, i)
        {
            const vector a = points[polygon[i]] - origin;
            const vector b = points[polygon[(i+1) % polygon.size()]] - origin;
            const scalar volume = orientation * (centre & (a ^ b)) / 6;
            total.volume += volume;
            total.first += 0.25 * volume * (centre + a + b);
        }
    }
}

bool identicalVertex(const point& a, const point& b)
{
    return a.x() == b.x() && a.y() == b.y() && a.z() == b.z();
}

bool samePolygon(const face& a, const face& b, const pointField& oldPoints,
    const pointField& newPoints, bool reversed = false)
{
    if (a.size() != b.size()) return false;
    forAll(b, start)
    {
        if (!identicalVertex(oldPoints[a[0]], newPoints[b[start]])) continue;
        bool equal = true;
        forAll(a, i)
        {
            const label index = reversed ? (start + b.size() - i) % b.size()
                : (start + i) % b.size();
            if (!identicalVertex(oldPoints[a[i]], newPoints[b[index]]))
            {
                equal = false;
                break;
            }
        }
        if (equal) return true;
    }
    return false;
}

bool overlaps(const point& low, const point& high, const point& a, const point& b)
{
    return low.x() < b.x() && high.x() > a.x()
        && low.y() < b.y() && high.y() > a.y()
        && low.z() < b.z() && high.z() > a.z();
}

bool logicalHex
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
            if (mag(t - 1) < allowedDisplacement) corner |= 1u << axis;
            else if (mag(t) >= allowedDisplacement) return false;
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
    boolList seen(oldFaces.size(), false);
    forAll(patches, patchi)
    {
        const polyPatch& patch = patches[patchi];
        if (patch.name() != names[patchi] || patch.type() != types[patchi]
            || patch.size() != sizes[patchi])
            FatalErrorInFunction << "Boundary patch changed: " << patch.name() << exit(FatalError);
        forAll(patch, i)
        {
            const label facei = patch.start() + i;
            const label source = mapping.faceMap()[facei] - boundaryStart;
            if (source < offset || source >= offset + sizes[patchi] || seen[source])
                FatalErrorInFunction << "Boundary-face map is not bijective within patch "
                    << patch.name() << ", new face " << facei << ", old face "
                    << mapping.faceMap()[facei] << exit(FatalError);
            seen[source] = true;
            if (!samePolygon(oldFaces[source], mesh.faces()[facei], oldPoints, mesh.points()))
                FatalErrorInFunction << "Boundary-face geometry changed on " << patch.name()
                    << ", new face " << facei << ", old face " << mapping.faceMap()[facei]
                    << ", vertex counts " << mesh.faces()[facei].size() << " and "
                    << oldFaces[source].size() << exit(FatalError);
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

    const IOdictionary protectedBox
    (
        IOobject("agglomerationBox", runTime.system(), mesh, IOobject::MUST_READ, IOobject::NO_WRITE)
    );
    const point protectedMin(protectedBox.get<point>("min"));
    const point protectedMax(protectedBox.get<point>("max"));

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
    std::map<Key, Octet, std::greater<Key>> groups;
    const labelListList& vertices = mesh.cellPoints();
    const scalarField oldVolumes(mesh.V());
    const vectorField& centres = mesh.cellCentres();
    const cellList& cells = mesh.cells();
    const faceList faces(mesh.faces());
    const labelList& owners = mesh.faceOwner();
    const labelList& neighbours = mesh.faceNeighbour();
    forAll(levels, celli)
    {
        const label level = levels[celli];
        if (level < 0 || level > 3 || cells[celli].size() != 6 || oldVolumes[celli] <= 0)
            continue;
        const scalar width = std::ldexp(20., -level);
        const vector grid = (centres[celli] - origin) / width;
        const label x = floor(grid.x()), y = floor(grid.y()), z = floor(grid.z());
        if (x < 0 || y < 0 || z < 0) continue;
        const point low = origin + width * vector(x, y, z);
        if (!logicalHex(vertices[celli], oldPoints, patchPoint, low, width)) continue;
        const point high = low + width * vector(1,1,1);
        // Protect the full support box plus at least two level-3 cells on each side.
        if (overlaps(low, high, protectedMin, protectedMax)) continue;
        const Key key{level, x/2, y/2, z/2};
        const label slot = x%2 + 2*(y%2) + 4*(z%2);
        Octet& octet = groups[key];
        if (octet.children[slot] != -1) octet.valid = false;
        octet.children[slot] = celli;
    }

    // Each logical pair differing in one axis must share exactly one quad.
    // Test projected levels with fine octets first. Keep distinct component
    // faces even when they end up connecting the same pair of cells.
    labelList masterForCell(oldCellCount, -1);
    labelList groupForCell(oldCellCount, -1);
    DynamicList<Parent> parents;
    labelList countByLevel(4, 0);
    for (const auto& entry : groups)
    {
        const Octet& octet = entry.second;
        const auto& children = octet.children;
        if (!octet.valid || std::find(children.begin(), children.end(), -1) != children.end())
            continue;
        const label childLevel = entry.first[0];
        bool valid = true;
        for (label slot = 0; slot < 8 && valid; ++slot)
        {
            const label child = children[slot];
            unsigned connections = 0;
            for (const label facei : cells[child])
            {
                if (facei >= boundaryStart || faces[facei].size() != 4)
                {
                    valid = false;
                    break;
                }
                const label other = owners[facei] == child ? neighbours[facei] : owners[facei];
                const auto found = std::find(children.begin(), children.end(), other);
                if (found != children.end())
                {
                    const label axisBit = slot ^ label(found - children.begin());
                    if ((axisBit != 1 && axisBit != 2 && axisBit != 4) || (connections & axisBit))
                    {
                        valid = false;
                        break;
                    }
                    connections |= axisBit;
                }
                else
                {
                    const label otherLevel = levels[other] - (masterForCell[other] >= 0);
                    if (mag(otherLevel - (childLevel - 1)) > 1)
                    {
                        valid = false;
                        break;
                    }
                }
            }
            valid = valid && connections == 7;
        }
        if (!valid) continue;
        ++countByLevel[childLevel];
        const label master = *std::min_element(children.begin(), children.end());
        Parent parent{0, point::zero, Moments()};
        for (const label child : children)
        {
            masterForCell[child] = master;
            groupForCell[child] = parents.size();
            parent.volume += oldVolumes[child];
            parent.centroid += oldVolumes[child] * centres[child];
        }
        parent.centroid /= parent.volume;
        for (const label child : children)
            addCellMoments(parent.exact, parent.centroid, cells[child], faces,
                oldPoints, owners, child);
        parents.append(parent);
    }
    const label nGroups = parents.size();
    Info << "Sibling octets by original level (0 through 3): " << countByLevel << nl;
    if (!nGroups)
        FatalErrorInFunction << "No eligible cell octets" << exit(FatalError);

    polyTopoChange topology(mesh);
    boolList usedPoint(mesh.nPoints(), false), retained(faces.size(), true);
    label removedFaces = 0;
    forAll(faces, facei)
    {
        const label owner = owners[facei];
        const label neighbour = facei < boundaryStart ? neighbours[facei] : -1;
        const label newOwner = masterForCell[owner] >= 0 ? masterForCell[owner] : owner;
        const label newNeighbour = neighbour >= 0 && masterForCell[neighbour] >= 0
            ? masterForCell[neighbour] : neighbour;
        if (newOwner == newNeighbour)
        {
            retained[facei] = false;
            topology.removeFace(facei, -1);
            ++removedFaces;
            continue;
        }
        for (const label vertex : faces[facei]) usedPoint[vertex] = true;
        if (newOwner == owner && newNeighbour == neighbour) continue;
        if (neighbour < 0)
            FatalErrorInFunction << "Selected octet touches a boundary face" << exit(FatalError);
        const label zone = mesh.faceZones().whichZone(facei);
        if (zone >= 0)
            FatalErrorInFunction << "Selected octet touches a face zone" << exit(FatalError);
        if (newOwner < newNeighbour)
            topology.modifyFace(faces[facei], facei, newOwner, newNeighbour, false, -1, -1, false);
        else
            topology.modifyFace(faces[facei].reverseFace(), facei, newNeighbour, newOwner, true,
                -1, -1, false);
    }
    if (removedFaces != 12*nGroups)
        FatalErrorInFunction << "Unexpected number of removed internal faces" << exit(FatalError);
    forAll(usedPoint, pointi)
        if (!usedPoint[pointi]) topology.removePoint(pointi, -1);
    forAll(masterForCell, celli)
        if (masterForCell[celli] >= 0 && masterForCell[celli] != celli)
            topology.removeCell(celli, masterForCell[celli]);

    autoPtr<mapPolyMesh> mapping = topology.changeMesh(mesh, false);
    mesh.updateMesh(mapping());
    if (mapping().hasMotionPoints()) mesh.movePoints(mapping().preMotionPoints());

    verifyPatches(mesh, mesh.boundaryMesh(), patchNames, patchTypes, patchSizes,
        oldFaces, oldPoints, boundaryStart, mapping());
    boolList seen(faces.size(), false);
    forAll(mesh.faces(), facei)
    {
        const label source = mapping().faceMap()[facei];
        if (source < 0 || source >= faces.size() || !retained[source] || seen[source])
            FatalErrorInFunction << "Non-bijective surviving face map at face " << facei
                << exit(FatalError);
        const bool flipped = mapping().flipFaceFlux().found(facei);
        if (!samePolygon(faces[source], mesh.faces()[facei], oldPoints, mesh.points(), flipped))
            FatalErrorInFunction << "Surviving face polygon changed at face " << facei
                << " (original " << source << ")" << exit(FatalError);
        seen[source] = true;
    }
    if (seen != retained)
        FatalErrorInFunction << "Surviving face set changed" << exit(FatalError);
    if (mesh.nCells() != oldCellCount - 7*nGroups)
        FatalErrorInFunction << "Unexpected resulting cell count" << exit(FatalError);

    const scalarField& volumes = mesh.V();
    labelList remapped(mesh.nCells());
    boolList visits(nGroups, false);
    scalar maxRelativeCentroidDifference = 0, maxAbsoluteCentroidDifference = 0;
    scalar maxExactVolumeDifference = 0, maxExactCentroidDifference = 0;
    label estimatorExceedances = 0;
    forAll(remapped, celli)
    {
        const label oldCell = mapping().cellMap()[celli];
        if (oldCell < 0 || oldCell >= oldCellCount)
            FatalErrorInFunction << "Result cell has no source cell" << exit(FatalError);
        const label group = groupForCell[oldCell];
        remapped[celli] = levels[oldCell] - (group >= 0);
        const scalar targetVolume = group >= 0 ? parents[group].volume : oldVolumes[oldCell];
        if (mag(volumes[celli] - targetVolume) > tolerance * max(1., mag(targetVolume)))
            FatalErrorInFunction << "Map disagrees with geometry at cell " << celli
                << ": expected volume " << targetVolume << ", got " << volumes[celli]
                << exit(FatalError);
        if (group >= 0)
        {
            const Parent& parent = parents[group];
            Moments exact;
            addCellMoments(exact, parent.centroid, mesh.cells()[celli],
                mesh.faces(), mesh.points(), mesh.faceOwner(), celli);
            const scalar side = std::cbrt(parent.exact.volume);
            const scalar exactVolumeDifference = mag(exact.volume - parent.exact.volume)
                / parent.exact.volume;
            const scalar exactCentroidDifference =
                mag(exact.first/exact.volume - parent.exact.first/parent.exact.volume) / side;
            maxExactVolumeDifference = max(maxExactVolumeDifference, exactVolumeDifference);
            maxExactCentroidDifference = max(maxExactCentroidDifference, exactCentroidDifference);
            if (visits[group] || exactVolumeDifference > 1e-9 || exactCentroidDifference > 1e-9)
                FatalErrorInFunction << "Integrated volume/centroid mismatch at cell " << celli
                    << ": relative volume " << exactVolumeDifference << ", relative centroid "
                    << exactCentroidDifference << exit(FatalError);

            // OpenFOAM's face-pyramid approximation can depend on its apex
            // for warped quads. Record it independently of the geometry gate.
            const scalar estimatorDifference = mag(mesh.cellCentres()[celli] - parent.centroid);
            maxAbsoluteCentroidDifference = max(maxAbsoluteCentroidDifference, estimatorDifference);
            maxRelativeCentroidDifference = max(maxRelativeCentroidDifference, estimatorDifference / side);
            estimatorExceedances += estimatorDifference / side > 1e-4;
            visits[group] = true;
        }
    }
    for (const bool visited : visits)
        if (!visited) FatalErrorInFunction << "Missing parent cell" << exit(FatalError);
    Info << "Maximum integrated relative volume difference: " << maxExactVolumeDifference << nl
         << "Maximum integrated centroid difference / parent side: " << maxExactCentroidDifference << nl
         << "Maximum OpenFOAM centroid estimator difference (m): " << maxAbsoluteCentroidDifference << nl
         << "Maximum OpenFOAM centroid estimator difference / parent side: "
         << maxRelativeCentroidDifference << nl
         << "OpenFOAM centroid estimator differences above 1e-4 of parent side: "
         << estimatorExceedances << nl;

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
