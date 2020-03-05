
#include "projet.hpp"

#ifdef SCENE_PROJECT

// Add vcl namespace within the current one - Allows to use function from vcl library without explicitely preceeding their name with vcl::
using namespace vcl;

mesh scene_model::create_water(int n_water)
{
    mesh water;

    for (int i = 0; i < n_water; i++)
    {
        for (int j = 0; j < n_water; j++)
        {

            water.position.push_back(vec3(i * 0.01f - 1, 1.0f, j * 0.01f - 1));
            water.color.push_back(vec4(0.5f, 0.5f, 1.0f, 1.0f));
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

            water.position.push_back(vec3(i * 0.01f - 1, 0.0f, j * 0.01f - 1));
            water.color.push_back(vec4(0.3f, 0.3f, 1.0f, 1.0f));
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

    return water;
}

mesh_drawable scene_model::create_glass(vec3 p0,vec3 u0,vec3 u1, vec3 u2) {

    mesh_drawable glass = mesh_drawable(mesh_primitive_parallelepiped(p0,u0,u1,u2));
    glass.uniform.color_alpha = 0.1f;
    glass.uniform.color = vec3(1.0f, 1.0f, 1.0f);
    glass.uniform.shading.specular = 0.8f;
    glass.uniform.shading.ambiant = 0.8f;
    glass.uniform.shading.diffuse = 0.8f;

    return glass;

}

void scene_model::setup_data(std::map<std::string, GLuint> &shaders, scene_structure &scene, gui_structure &gui)
{

    // Initial position of the camera
    scene.camera.translation = {0.0f, -0.5f, 0.0f};
    scene.camera.scale = 5.0f;
    gui.show_frame_camera = false;

    //Water

    int n_water = 200;
    water = create_water(n_water);

    water_draw = water;
    water_draw.shader = shaders["mesh"];
    water_draw.uniform.color_alpha = 0.3f;
    water_draw.uniform.shading.specular = 0.6f;
    water_draw.uniform.shading.ambiant = 0.7f;
    water_draw.uniform.shading.diffuse = 0.5f;

    //Glass
    glass_back = create_glass({-1.05, -0.05, -1.05}, {2.1, 0, 0}, {0, 1.25, 0}, {0, 0, 0.05});
    glass_front = create_glass({-1.05, -0.05, 1}, {2.1, 0, 0}, {0, 1.25, 0}, {0, 0, 0.05});
    glass_left = create_glass({-1.05, -0.05, -1.05}, {0.05, 0, 0}, {0, 1.25, 0}, {0, 0, 2.1});
    glass_right = create_glass({1, -0.05, -1.05}, {0.05, 0, 0}, {0, 1.25, 0}, {0, 0, 2.1});
    glass_bottom = create_glass({-1.05, -0.05, -1.05}, {2.1, 0, 0}, {0, 0.05, 0}, {0, 0, 2.1});

    cube = mesh_drawable(mesh_primitive_parallelepiped());

    cube.shader = shaders["mesh"];
    cube.uniform.color_alpha = 1.0f;
    cube.uniform.color = vec3(1.0f, 0.0f, 0.0f);
    cube.uniform.shading.specular = 0.8f;
    cube.uniform.shading.ambiant = 0.5f;
    cube.uniform.shading.diffuse = 0.8f;
    cube.uniform.transform.scaling = 0.3f;
    cube.uniform.transform.translation = {0, 0.5, 0.5};

    sphere = mesh_drawable(mesh_primitive_sphere(1.0f));
	sphere.shader = shaders["mesh"];
	init_fishes();

}

void scene_model::update_water(float dt)
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

    if (c)
    {
        draw(cube, scene.camera);
        if (gui_scene.wireframe)
            draw(cube, scene.camera, shaders["wireframe"]);
    }

    if (f) {

    display_fishes(scene);

    }

    set_gui();
    update_water(dt);
    compute_time_step_fishes(dt);

    glEnable(GL_BLEND);
    glDepthMask(false);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);


    if (g)
    {
        draw(glass_back, scene.camera, shaders["mesh"]);
        if (gui_scene.wireframe)
            draw(glass_back, scene.camera, shaders["wireframe"]);

        draw(glass_front, scene.camera, shaders["mesh"]);
        if (gui_scene.wireframe)
            draw(glass_front, scene.camera, shaders["wireframe"]);

        draw(glass_left, scene.camera, shaders["mesh"]);
        if (gui_scene.wireframe)
            draw(glass_left, scene.camera, shaders["wireframe"]);

        draw(glass_right, scene.camera, shaders["mesh"]);
        if (gui_scene.wireframe)
            draw(glass_right, scene.camera, shaders["wireframe"]);

        draw(glass_bottom, scene.camera, shaders["mesh"]);
        if (gui_scene.wireframe)
            draw(glass_bottom, scene.camera, shaders["wireframe"]);
    }

    if (w)
    {
        draw(water_draw, scene.camera);
        if (gui_scene.wireframe)
            draw(water_draw, scene.camera, shaders["wireframe"]);
    }

    glDepthMask(true);
}

void scene_model::set_gui()
{
    ImGui::Checkbox("Wireframe", &gui_scene.wireframe);
    ImGui::Checkbox("Water", &w);
    ImGui::Checkbox("Glass", &g);
    ImGui::Checkbox("Cube", &c);
    ImGui::Checkbox("Fish", &f);
}


void scene_model::init_fishes() {
	for (size_t k = 0; k < N_fishes; k++) {
		particle_structure fish;

		// Initial position
		float x = rand_interval(0, 0.5);
		float y = rand_interval(0, 0.5);
		float z = rand_interval(0, 0.5);
		fish.p = vec3(x, y, z);
		std::cout << fish.p << std::endl;

		// Initial speed
		float vx = rand_interval(0, 1);
		float vy = rand_interval(0, 1);
		float vz = rand_interval(0, 1);
		fish.v = vec3(vx, vy, vz);
		fish.f = vec3(0,0,0);

		fishes.push_back(fish);
	}


}
void scene_model::display_fishes(scene_structure& scene)
{
	const size_t N = fishes.size();
	for (size_t k = 0; k < N; ++k)
	{
		const particle_structure& part = fishes[k];

		sphere.uniform.transform.translation = part.p;
		sphere.uniform.transform.scaling = 0.05f;
		vec3 c = { 1,0,0 };
		sphere.uniform.color = c;
		draw(sphere, scene.camera);
	}
}

float f(float x) {
	//return 0 if too long
	float r = 3.0;
	//if (x > r) {
		//return 0;
	//}
	//float alpha1 = 0.1;
	//float alpha2 = 0.000001;
	//return 0.001* alpha1 / (x * x) - alpha2 / (x * x * x * x);
	float alpha1 = 0.08;
	float alpha2 = 0.0005;
	float x0 = 0.1;
	float k = 0.01;
	return alpha1 * exp(-k * (pow((x - x0) / x0, 2 ))) - alpha2 * exp(- x / x0);

}
void scene_model::compute_time_step_fishes(float dt)
{
	// Set forces
	const size_t N = fishes.size();
	/*
	for (size_t k = 0; k < N; ++k) {
		//frottement fluide
		float alpha = 0.01;
		fishes[k].f -= alpha * fishes[k].v;

		//boids
		for (size_t j = 0; j < N; j++) {

			float d = norm(fishes[j].p - fishes[k].p);
			if (d > 1e-5) {
				fishes[k].f += f(d) * (fishes[j].p - fishes[k].p) / d;
			}
		}


	}
	*/
	float rsep = 0.5;
	float ralign = 1.0;
	float rcohez = 2.0;
	float vmax = 5.0f;
	for (size_t k = 0; k < N; ++k) {
		vec3 vsep = { 0,0,0 };
		vec3 valign = { 0,0,0 };
		int nalign = 0;
		int ncohez = 0;
		vec3 vcohez = { 0,0,0 };
		particle_structure& fish = fishes[k];
		vec3& v = fish.v;
		vec3& p = fish.p;

		//frottement
		float alpha_frot = 0.3;
		fish.f = -alpha_frot*fish.v;

		//attraction to the center
        vec3 center = {0,0.5,0};
		float raideur = 0.1;
		float l0 = 0.1f;
		fish.f += - raideur * (norm(p - center)-l0) * (fish.p - center);

		//boids model
		for (size_t j = 0; j < N; j++) {
			particle_structure neighbor = fishes[j];
			float dist = norm(fish.p - neighbor.p);
			if (dist > 0) {
				//separation
				if (dist < rsep) {
					vsep += (dist / rsep - 1) / dist * (neighbor.p - fish.p);
				}
				//alignment
				if (dist < ralign) {
					valign += neighbor.v;
					nalign += 1;
				}
				//cohesion
				if (dist < rcohez) {
					vcohez += (neighbor.p - fish.p) / dist;
				}
			}
			
		}
		if (nalign > 0) {
			valign = valign / nalign;
		}
		if (ncohez > 0) {
			vcohez = vcohez / ncohez;
		}

		float alpha_sep = 0.1;
		float alpha_align = 0.02;
		float alpha_cohez = 0.02;
		v += alpha_sep * vsep + alpha_align * valign + alpha_align * vcohez;
		if (norm(v) > vmax) {
			v = v * vmax / norm(v);
		}
		//perturbation al�atoire
		float alpha_alea = 0.1;
		float vx = rand_interval(-1, 1);
		float vy = rand_interval(-1, 1);
		float vz = rand_interval(-1, 1);
		fish.v += alpha_alea * vec3(vx, vy, vz);

		//force repulsive bords :
		float drepuls = 0.1f;
		float frepuls = 50.0f;
		if (p.x > 1-drepuls) {
			fish.f += vec3(-frepuls, 0, 0);
		}
		if (p.x < -1+drepuls) {
			fish.f += vec3(frepuls, 0, 0);
		}
		if (p.y > 1 - drepuls) {
			fish.f += vec3(0, -frepuls, 0);
		}
		if (p.y < drepuls) {
			fish.f += vec3(0, frepuls, 0);
		}
		if (p.z > 1 - drepuls) {
			fish.f += vec3(0, 0, -frepuls);
		}
		if (p.z < -1+drepuls) {
			fish.f += vec3(0, 0, frepuls);
		}
	}

	// Integrate position and speed of particles through time
	for (size_t k = 0; k < N; ++k) {
		particle_structure& fish = fishes[k];
		vec3& v = fish.v;
		vec3& p =fish.p;
		vec3 const& f = fish.f;

		v = v + dt * f; 
		p = p + dt * v;
	}
}

#endif
