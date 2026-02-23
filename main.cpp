#include "geogram/mesh/mesh_io.h"
#include "untangler3d.hpp"
#include <geogram/mesh/mesh.h>


using namespace GEO;

int main(int argc, char *argv[])
{
    GEO::initialize();

    std::string reference_file = argv[1];
    std::string deformed_file = argv[2];

    std::cout << reference_file << std::endl;
    std::cout << deformed_file << std::endl;

    Mesh reference_mesh;
    if (!mesh_load(reference_file, reference_mesh))
    {
        std::cerr << "Error loading mesh: " << deformed_file << std::endl;
        exit(-1);
    }

    Mesh deformed_mesh;
    if (!mesh_load(deformed_file, deformed_mesh))
    {
        std::cerr << "Error loading mesh: " << deformed_file << std::endl;
        exit(-1);
    }

    Untangler3D untangler(reference_mesh, deformed_mesh);

    {
        Timer t("TIMER go: ");
        untangler.go();
    }

    MeshIOFlags flags;
    OutputGeoFile geofile("out.geogram");
    mesh_save(deformed_mesh, geofile, flags);


    return 0;
}
