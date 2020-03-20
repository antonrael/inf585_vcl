#pragma once

#include "scenes/base/base.hpp"

#ifdef SCENE_PROJECT


struct gui_scene_structure
{
	bool wireframe = false; // Display the wireframe
	int face_index = 0;

};

struct particle_structure
{
	vcl::vec3 p; // Position
	vcl::vec3 v; // Speed
	vcl::vec3 f; // Forces

};


struct scene_model : scene_base
{

	vcl::mesh create_water(int n_water);
	vcl::mesh create_sand(int n_sand);

	int nb_seaweeds = 200;
	int size_seaweed = 10;
	vcl::buffer<vcl::mesh> create_seaweeds(int nb_seaweeds, std::uint32_t size_seaweed);


	vcl::mesh_drawable create_glass(vcl::vec3 p0, vcl::vec3 u0, vcl::vec3 u1, vcl::vec3 u2);
	void setup_data(std::map<std::string, GLuint>& shaders, scene_structure& scene, gui_structure& gui);
	void frame_draw(std::map<std::string, GLuint>& shaders, scene_structure& scene, gui_structure& gui);

	void set_gui();

	void update_water(float dt);
	void update_seaweeds(int nb_seaweeds, float dt);

	void compute_forces_seaweeds(int number_of_seaweeds, int n_seaweed);
	void update_seaweeds(int nb_seaweeds, int size_seaweed, float dt);

	vcl::mesh water;
	vcl::mesh_drawable water_draw;

	vcl::mesh sand;
	vcl::mesh_drawable sand_draw;


	vcl::mesh_drawable glass_left;
	vcl::mesh_drawable glass_right;
	vcl::mesh_drawable glass_front;
	vcl::mesh_drawable glass_back;
	vcl::mesh_drawable glass_bottom;

	vcl::buffer<vcl::mesh> seaweeds;
	vcl::buffer<vcl::buffer<vcl::vec3>> initial_positions_seaweeds;
	vcl::buffer<vcl::buffer<vcl::vec3>> speeds_seaweeds;
	vcl::buffer<vcl::buffer<vcl::vec3>> forces_seaweeds;
	vcl::buffer<vcl::mesh_drawable> seaweeds_draw;

	vcl::mesh_drawable cube;


	//Fishes
	int N_fishes = 40;
	void init_fishes();
	void compute_time_step_fishes(float dt);
	std::vector<particle_structure> fishes;
	vcl::mesh_drawable sphere;      // Visual display of particles
	//display fish hierarchy
	std::vector<vcl::hierarchy_mesh_drawable> hierarchy_fishes;
	vcl::hierarchy_mesh_drawable create_fish(std::map<std::string, GLuint>& shaders);
	void display_fish(std::map<std::string, GLuint>& shaders, scene_structure& scene);

	//textures

	GLuint texture_sand;
	GLuint texture_glass;


	vcl::timer_event timer;
	gui_scene_structure gui_scene;

	bool w = true;
	bool c = false;
	bool g = true;
	bool f = true;
	bool s = true;
	bool sw = true;

};

#endif