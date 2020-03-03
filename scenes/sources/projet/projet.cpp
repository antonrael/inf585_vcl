
#include "projet.hpp"

#ifdef SCENE_PROJECT

// Add vcl namespace within the current one - Allows to use function from vcl library without explicitely preceeding their name with vcl::
using namespace vcl;

void scene_model::setup_data(std::map<std::string, GLuint> &shaders, scene_structure &scene, gui_structure &gui)
{

    int n_water = 200;

    // Initial position of the camera
    scene.camera.translation = {0.0f, -0.5f, 0.0f};
    scene.camera.scale = 5.0f;
    gui.show_frame_camera = false;

    for (int i = 0; i < n_water; i++)
    {
        for (int j = 0; j < n_water; j++)
        {

            water.position.push_back(vec3(i * 0.01f - 1, 1.0f, j * 0.005f));
            water.color.push_back(vec4(0.5f, 0.5f, 1.0f, 0.0f));
        }
    }

    for (uint i = 0; i < n_water; i++)
    {
        for (uint j = 0; j < n_water; j++)
        {

            if (i > 0 && j > 0)
            {
                water.connectivity.push_back({i * n_water + j, (i - 1) * n_water + j, i * n_water + j - 1});
            }

            if (i < n_water - 1 && j < n_water - 1)
            {
                water.connectivity.push_back({i * n_water + j, (i + 1) * n_water + j, i * n_water + j + 1});
            }
        }
    }

    for (int i = 0; i < n_water; i++)
    {
        for (int j = 0; j < n_water; j++)
        {

            water.position.push_back(vec3(i * 0.01f - 1, 0.0f, j * 0.005f));
            water.color.push_back(vec4(0.3f, 0.3f, 1.0f, 0.0f));
        }
    }

    for (uint i = 0; i < n_water; i++)
    {
        for (uint j = 0; j < n_water; j++)
        {

            if (i > 0 && j > 0)
            {
                water.connectivity.push_back({n_water * n_water + i * n_water + j, n_water * n_water + (i - 1) * n_water + j, n_water * n_water + i * n_water + j - 1});
            }

            if (i < n_water - 1 && j < n_water - 1)
            {
                water.connectivity.push_back({n_water * n_water + i * n_water + j, n_water * n_water + (i + 1) * n_water + j, n_water * n_water + i * n_water + j + 1});
            }
        }
    }

    for (uint j = 0; j < n_water - 1; j++)
    {
        water.connectivity.push_back({j, j + 1, n_water * n_water + j});
        water.connectivity.push_back({j + n_water * (n_water - 1), j + 1 + n_water * (n_water - 1), n_water * n_water + j + n_water * (n_water - 1)});
        water.connectivity.push_back({j + 1, n_water * n_water + j + 1, n_water * n_water + j});
        water.connectivity.push_back({j + 1 + n_water * (n_water - 1), n_water * n_water + j + 1 + n_water * (n_water - 1), n_water * n_water + j + n_water * (n_water - 1)});
    }

    for (uint i = 0; i < n_water - 1; i++)
    {
        water.connectivity.push_back({i * n_water, (i + 1) * n_water, n_water * n_water + i * n_water});
        water.connectivity.push_back({i * n_water + n_water - 1, (i + 1) * n_water + n_water - 1, n_water * n_water + i * n_water + n_water - 1});

        water.connectivity.push_back({(i + 1) * n_water, (i + 1) * n_water + n_water * n_water, n_water * n_water + i * n_water});
        water.connectivity.push_back({(i + 1) * n_water + n_water - 1, (i + 1) * n_water + n_water * n_water + n_water - 1, n_water * n_water + i * n_water + n_water - 1});
    }

    normal(water.position, water.connectivity, water.normal, true);

    water_draw = water;
    water_draw.shader = shaders["mesh"];
    water_draw.uniform.color_alpha = 1.0f;
    water_draw.uniform.shading.specular = 0.3f;
    water_draw.uniform.shading.ambiant = 0.7f;
    water_draw.uniform.shading.diffuse = 0.5f;

    //Outside glass

    glass.position.push_back(vec3(-1.05f, -0.05f, -0.05f));
    glass.position.push_back(vec3(-1.05f, -0.05f, 1.05f));
    glass.position.push_back(vec3(1.05f, -0.05f, -0.05f));
    glass.position.push_back(vec3(1.05f, -0.05f, 1.05f));

    glass.position.push_back(vec3(-1.05f, 1.2f, -0.05f));
    glass.position.push_back(vec3(-1.05f, 1.2f, 1.05f));
    glass.position.push_back(vec3(1.05f, 1.2f, -0.05f));
    glass.position.push_back(vec3(1.05f, 1.2f, 1.05f));

    //Inside glass

    glass.position.push_back(vec3(-1.0f, 0.0f, 0.0f));
    glass.position.push_back(vec3(-1.0f, 0.0f, 1.0f));
    glass.position.push_back(vec3(1.0f, 0.0f, 0.0f));
    glass.position.push_back(vec3(1.0f, -0.0f, 1.0f));

    glass.position.push_back(vec3(-1.0f, 1.2f, 0.0f));
    glass.position.push_back(vec3(-1.0f, 1.2f, 1.0f));
    glass.position.push_back(vec3(1.0f, 1.2f, 0.0f));
    glass.position.push_back(vec3(1.0f, 1.2f, 1.0f));

    //Connect Outside Glass

    glass.connectivity.push_back({0, 1, 2});
    glass.connectivity.push_back({1, 2, 3});
    glass.connectivity.push_back({0, 1, 4});
    glass.connectivity.push_back({1, 4, 5});
    glass.connectivity.push_back({0, 2, 4});
    glass.connectivity.push_back({2, 4, 6});
    glass.connectivity.push_back({1, 3, 5});
    glass.connectivity.push_back({3, 5, 7});
    glass.connectivity.push_back({2, 3, 6});
    glass.connectivity.push_back({3, 6, 7});

    //Connect Inside Glass

    glass.connectivity.push_back({8, 9, 10});
    glass.connectivity.push_back({9, 10, 11});
    glass.connectivity.push_back({8, 9, 12});
    glass.connectivity.push_back({9, 12, 13});
    glass.connectivity.push_back({8, 10, 12});
    glass.connectivity.push_back({10, 12, 14});
    glass.connectivity.push_back({9, 11, 13});
    glass.connectivity.push_back({11, 13, 15});
    glass.connectivity.push_back({10, 11, 14});
    glass.connectivity.push_back({11, 14, 15});

    //Connect Inside and Oustide Glass

    glass.connectivity.push_back({4, 12, 5});
    glass.connectivity.push_back({5, 12, 13});
    glass.connectivity.push_back({5, 13, 15});
    glass.connectivity.push_back({5, 15, 7});
    glass.connectivity.push_back({4, 6, 12});
    glass.connectivity.push_back({6, 12, 14});
    glass.connectivity.push_back({6, 14, 15});
    glass.connectivity.push_back({6, 7, 15});

    normal(glass.position, glass.connectivity, glass.normal);

    glass_draw = glass;
    glass_draw.shader = shaders["mesh"];
    glass_draw.uniform.color_alpha = 0.0f;
    glass_draw.uniform.color = vec3(0.8f, 0.9f, 0.9f);
    glass_draw.uniform.shading.specular = 0.8f;
    glass_draw.uniform.shading.ambiant = 0.8f;
    glass_draw.uniform.shading.diffuse = 0.0f;
}

void scene_model::compute_time_step(float dt)
{

    int n_water = std::sqrt(water.position.size() / 2);

    for (size_t i = 0; i < n_water; ++i)
    {

        for (size_t j = 0; j < n_water; ++j)
        {
            int k = i * n_water + j;
            float x = 0.02 * i;
            float y = 0.02 * j;
            float t = timer.t;
            water.position[k][1] = 1.0f + 0.07 * perlin(0.5 * x, 0.5 * y, t, 2, 0.3);
        }
    }

    normal(water.position, water.connectivity, water.normal, true);

    water_draw.update_position(water.position);
    water_draw.update_normal(water.normal);
}

void scene_model::frame_draw(std::map<std::string, GLuint> &shaders, scene_structure &scene, gui_structure &)
{

    float dt = 0.02f;
    if (!(timer.update() > 0))
        dt = 0;

    set_gui();
    compute_time_step(dt);

    glEnable(GL_BLEND);

    //glDepthMask(false);

    glBlendFunc(GL_DST_COLOR, GL_ONE_MINUS_DST_COLOR);

    draw(glass_draw, scene.camera);
    if (gui_scene.wireframe)
        draw(glass_draw, scene.camera, shaders["wireframe"]);

    glBlendFunc(GL_DST_COLOR, GL_ONE_MINUS_DST_COLOR);

    draw(water_draw, scene.camera);
    if (gui_scene.wireframe)
        draw(water_draw, scene.camera, shaders["wireframe"]);
}

void scene_model::set_gui()
{
    ImGui::Checkbox("Wireframe", &gui_scene.wireframe);
}

#endif
