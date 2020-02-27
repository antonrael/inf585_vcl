#pragma once

#include "scenes/base/base.hpp"

#ifdef SCENE_PROJECT


struct gui_scene_structure
{
    bool wireframe   = false; // Display the wireframe
    int face_index = 0;

};



struct scene_model : scene_base
{

    void setup_data(std::map<std::string,GLuint>& shaders, scene_structure& scene, gui_structure& gui);
    void frame_draw(std::map<std::string,GLuint>& shaders, scene_structure& scene, gui_structure& gui);

    void set_gui();

    void compute_time_step(float dt);

    vcl::mesh water;
    vcl::mesh_drawable water_draw;


    vcl::timer_event timer;
    gui_scene_structure gui_scene;

};

#endif