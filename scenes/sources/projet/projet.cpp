
#include "projet.hpp"

#ifdef SCENE_PROJECT

// Add vcl namespace within the current one - Allows to use function from vcl library without explicitely preceeding their name with vcl::
using namespace vcl;

void scene_model::setup_data(std::map<std::string, GLuint> &shaders, scene_structure &scene, gui_structure &gui)
{

    int n_water = 200;

    // Initial position of the camera
    scene.camera.translation = {0.0f, 0.0f, 0.0f};
    scene.camera.scale = 5.0f;
    gui.show_frame_camera = false;

    for (int i = 0; i < n_water; i++)
    {
        for (int j = 0; j < n_water; j++)
        {

            water.position.push_back(vec3(i * 0.005f, 1.0f, j * 0.005f));
            water.color.push_back(vec4(0.6f, 0.6f, 1.0f, 0.0f));
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

            water.position.push_back(vec3(i * 0.005f, 0.0f, j * 0.005f));
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

    normal(water.position, water.connectivity, water.normal);

    water_draw = water;
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
            //std::cout << x << " " << y << std::endl;
            float t = timer.t;
            water.position[k][1] = 1.0f + 0.05 * perlin(x, y, t, 2, 0.3);
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

    draw(water_draw, scene.camera, shaders["water"]);
}

void scene_model::set_gui()
{
    ImGui::Checkbox("Wireframe", &gui_scene.wireframe);
}

#endif
