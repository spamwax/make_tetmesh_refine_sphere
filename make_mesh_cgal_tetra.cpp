#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
//#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Mesh_3/Robust_intersection_traits_3.h>

#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>

#include <CGAL/Polyhedral_mesh_domain_3.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/refine_mesh_3.h>

#include <CGAL/IO/Polyhedron_iostream.h>
#include <fstream>

// Computational Kernel Definitions necessary for robust meshing
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
//typedef CGAL::Exact_predicates_exact_constructions_kernel K;
typedef CGAL::Mesh_3::Robust_intersection_traits_3<K> Geom_traits;
typedef CGAL::Polyhedron_3<Geom_traits> Polyhedron;
//typedef CGAL::Polyhedron_3<K>      Polyhedron;
typedef CGAL::Polyhedral_mesh_domain_3<Polyhedron, Geom_traits> Mesh_domain;
//typedef CGAL::Polyhedral_mesh_domain_3<Polyhedron> Mesh_domain;

// Meshing typedefs
typedef CGAL::Mesh_triangulation_3<Mesh_domain>::type Tr;
typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;

// Parameters
typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;
using namespace CGAL::parameters;

int main(int argc, char *argv[]) {
	Polyhedron cutbrain;
	std::string ifn, cfn;
	double facet_angle=25, facet_size=3.0, facet_distance=0.5,
	       cell_radius_edge=3,cell_size=3.0;
	bool defulatcriteria = false;

	if (argc == 1) {
		std::cout << " Enter the surface mesh filename: ";
		std:: cin >> ifn;
		defulatcriteria = true;
		std::cout << " (Using default settings for meshing parameters!)\n";
	}
	else if (argc==2) {
		ifn = argv[1];
		defulatcriteria = true;
	}
	else if (argc == 3) {
		ifn = argv[1];
		cfn = argv[2];
	}

	std::ifstream ifs(ifn.c_str());
	if (!ifs) {
		std::cerr << " Could not open input surface mesh file!" << std::endl;
		exit(EXIT_FAILURE);
	}

	if (!defulatcriteria) {
        std::cout << " Reading criteria from file." << std::endl;
		std::ifstream cfs(cfn.c_str());
		if (!cfs) {
			std::cerr << " Can not read mesh criteria file!\n";
			exit(-1);
		}
		cfs >> facet_angle;
		cfs >> facet_size;
		cfs >> facet_distance;
		cfs >> cell_radius_edge;
		cfs >> cell_size;
	}

	std::cout << "    cgal_mesh: reading the input polyhedron.\n";
	ifs >> cutbrain;
	std::cout << "    cgal_mesh: done reading the input polyhedron.\n";
	
	if (!cutbrain.is_closed()) {
		std::cerr << " Input polyhedron is not closed!\n";
		exit(-1);
	}
	
	Mesh_domain domain(cutbrain);

	std::cout << "    cgal_mesh: setting domain criteria.\n";
	Mesh_criteria criteria(facet_angle, facet_size, facet_distance, cell_radius_edge, cell_size);

	std::cout << "    cgal_mesh: running make_mesh_3.\n";
	C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria);

	std::ofstream medit_file("out.mesh");
	c3t3.output_to_medit(medit_file);
	medit_file.close();
	std::cout << facet_angle << ' ' << std::endl;
	std::cout << facet_size << ' ' << std::endl;
	std::cout << facet_distance << ' ' << std::endl;
	std::cout << cell_radius_edge << ' ' << std::endl;
	std::cout << cell_size << ' ' << std::endl;
	return EXIT_SUCCESS;
	
}

