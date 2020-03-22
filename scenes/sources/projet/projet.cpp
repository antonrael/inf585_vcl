
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

	for (std::uint32_t i = 0; i < n_water; i++)
	{
		for (std::uint32_t j = 0; j < n_water; j++)
		{

			if (i > 0 && j > 0)
			{
				water.connectivity.push_back({ i * n_water + j, (i - 1) * n_water + j, i * n_water + j - 1 });
			}

			if (i < n_water - 1 && j < n_water - 1)
			{
				water.connectivity.push_back({ i * n_water + j, (i + 1) * n_water + j, i * n_water + j + 1 });
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

	for (std::uint32_t i = 0; i < n_water; i++)
	{
		for (std::uint32_t j = 0; j < n_water; j++)
		{

			if (i > 0 && j > 0)
			{
				water.connectivity.push_back({ n_water * n_water + i * n_water + j, n_water * n_water + (i - 1) * n_water + j, n_water * n_water + i * n_water + j - 1 });
			}

			if (i < n_water - 1 && j < n_water - 1)
			{
				water.connectivity.push_back({ n_water * n_water + i * n_water + j, n_water * n_water + (i + 1) * n_water + j, n_water * n_water + i * n_water + j + 1 });
			}
		}
	}

	for (std::uint32_t j = 0; j < n_water - 1; j++)
	{
		water.connectivity.push_back({ j, j + 1, n_water * n_water + j });
		water.connectivity.push_back({ j + n_water * (n_water - 1), j + 1 + n_water * (n_water - 1), n_water * n_water + j + n_water * (n_water - 1) });
		water.connectivity.push_back({ j + 1, n_water * n_water + j + 1, n_water * n_water + j });
		water.connectivity.push_back({ j + 1 + n_water * (n_water - 1), n_water * n_water + j + 1 + n_water * (n_water - 1), n_water * n_water + j + n_water * (n_water - 1) });
	}

	for (std::uint32_t i = 0; i < n_water - 1; i++)
	{
		water.connectivity.push_back({ i * n_water, (i + 1) * n_water, n_water * n_water + i * n_water });
		water.connectivity.push_back({ i * n_water + n_water - 1, (i + 1) * n_water + n_water - 1, n_water * n_water + i * n_water + n_water - 1 });

		water.connectivity.push_back({ (i + 1) * n_water, (i + 1) * n_water + n_water * n_water, n_water * n_water + i * n_water });
		water.connectivity.push_back({ (i + 1) * n_water + n_water - 1, (i + 1) * n_water + n_water * n_water + n_water - 1, n_water * n_water + i * n_water + n_water - 1 });
	}

	normal(water.position, water.connectivity, water.normal, true);

	return water;
}

mesh scene_model::create_sand(int n_sand)
{

	mesh sand;

	for (int i = 0; i < n_sand; i++)
	{
		for (int j = 0; j < n_sand; j++)
		{
			sand.position.push_back(vec3(i * 0.01f - 1, 0.05f + 0.02 * perlin(0.1 * i, 0.1 * j, 4, 0.3), j * 0.01f - 1));
			sand.color.push_back(vec4(1.0f, 1.0f, 1.0f, 1.0f));
			sand.texture_uv.push_back({ float(i * 0.02), float(j * 0.02) });
		}
	}

	for (std::uint32_t i = 0; i < n_sand; i++)
	{
		for (std::uint32_t j = 0; j < n_sand; j++)
		{

			if (i > 0 && j > 0)
			{
				sand.connectivity.push_back({ i * n_sand + j, (i - 1) * n_sand + j, i * n_sand + j - 1 });
			}

			if (i < n_sand - 1 && j < n_sand - 1)
			{
				sand.connectivity.push_back({ i * n_sand + j, (i + 1) * n_sand + j, i * n_sand + j + 1 });
			}
		}
	}

	for (int i = 0; i < n_sand; i++)
	{
		for (int j = 0; j < n_sand; j++)
		{

			sand.position.push_back(vec3(i * 0.01f - 1, 0.005f, j * 0.01f - 1));
			sand.color.push_back(vec4(0.3f, 0.3f, 1.0f, 1.0f));
			sand.texture_uv.push_back({ float((i + 5) * 0.02), float((j - 5) * 0.02) });
		}
	}

	for (std::uint32_t i = 0; i < n_sand; i++)
	{
		for (std::uint32_t j = 0; j < n_sand; j++)
		{

			if (i > 0 && j > 0)
			{
				sand.connectivity.push_back({ n_sand * n_sand + i * n_sand + j, n_sand * n_sand + (i - 1) * n_sand + j, n_sand * n_sand + i * n_sand + j - 1 });
			}

			if (i < n_sand - 1 && j < n_sand - 1)
			{
				sand.connectivity.push_back({ n_sand * n_sand + i * n_sand + j, n_sand * n_sand + (i + 1) * n_sand + j, n_sand * n_sand + i * n_sand + j + 1 });
			}
		}
	}

	for (std::uint32_t j = 0; j < n_sand - 1; j++)
	{
		sand.connectivity.push_back({ j, j + 1, n_sand * n_sand + j });
		sand.connectivity.push_back({ j + n_sand * (n_sand - 1), j + 1 + n_sand * (n_sand - 1), n_sand * n_sand + j + n_sand * (n_sand - 1) });
		sand.connectivity.push_back({ j + 1, n_sand * n_sand + j + 1, n_sand * n_sand + j });
		sand.connectivity.push_back({ j + 1 + n_sand * (n_sand - 1), n_sand * n_sand + j + 1 + n_sand * (n_sand - 1), n_sand * n_sand + j + n_sand * (n_sand - 1) });
	}

	for (std::uint32_t i = 0; i < n_sand - 1; i++)
	{
		sand.connectivity.push_back({ i * n_sand, (i + 1) * n_sand, n_sand * n_sand + i * n_sand });
		sand.connectivity.push_back({ i * n_sand + n_sand - 1, (i + 1) * n_sand + n_sand - 1, n_sand * n_sand + i * n_sand + n_sand - 1 });

		sand.connectivity.push_back({ (i + 1) * n_sand, (i + 1) * n_sand + n_sand * n_sand, n_sand * n_sand + i * n_sand });
		sand.connectivity.push_back({ (i + 1) * n_sand + n_sand - 1, (i + 1) * n_sand + n_sand * n_sand + n_sand - 1, n_sand * n_sand + i * n_sand + n_sand - 1 });
	}

	normal(sand.position, sand.connectivity, sand.normal, true);

	return sand;
}

buffer<vcl::mesh> scene_model::create_seaweeds(int nb_seaweeds, std::uint32_t size_seaweed)
{

	buffer<vcl::mesh> seaweeds;

	std::default_random_engine generator;
	std::uniform_real_distribution<double> distribution(-0.9, 0.9);

	double height = 0.5f;
	double width = 0.006f;

	for (size_t i = 0; i < nb_seaweeds; i++)
	{

		double delta_height = 0.1 * distribution(generator);
		float posx = distribution(generator);
		float posz = distribution(generator);
		double delta_width = 0.00 * distribution(generator);

		mesh shape;

		for (int i = 0; i < size_seaweed + 1; i++)
		{

			shape.position.push_back({ posx, static_cast<float>(((float)i / (float)size_seaweed) * (float)(height + delta_height) + 0.04f), posz });
			shape.position.push_back({ posx, static_cast<float>(((float)i / (float)size_seaweed) * (float)(height + delta_height) + 0.04f), posz + static_cast<float>(width + delta_width) });
			shape.position.push_back({ posx + static_cast<float>(width + delta_width), static_cast<float>(((float)i / (float)size_seaweed) * (float)(height + delta_height) + 0.04f), posz + static_cast<float>(width + delta_width) });
			shape.position.push_back({ posx + static_cast<float>(width + delta_width), static_cast<float>(((float)i / (float)size_seaweed) * (float)(height + delta_height) + 0.04f), posz });
		}

		for (std::uint32_t i = 0; i < size_seaweed; i++)
		{

			shape.connectivity.push_back({ 4 * i, 4 * i + 1, 4 * i + 5 });
			shape.connectivity.push_back({ 4 * i, 4 * i + 4, 4 * i + 5 });
			shape.connectivity.push_back({ 4 * i + 1, 4 * i + 2, 4 * i + 5 });
			shape.connectivity.push_back({ 4 * i + 6, 4 * i + 2, 4 * i + 5 });
			shape.connectivity.push_back({ 4 * i + 3, 4 * i + 2, 4 * i + 7 });
			shape.connectivity.push_back({ 4 * i + 6, 4 * i + 2, 4 * i + 7 });
			shape.connectivity.push_back({ 4 * i + 3, 4 * i, 4 * i + 7 });
			shape.connectivity.push_back({ 4 * i + 4, 4 * i, 4 * i + 7 });
		}

		shape.connectivity.push_back({ 0, 1, 2 });
		shape.connectivity.push_back({ 3, 1, 2 });
		shape.connectivity.push_back({ 4 * size_seaweed, 4 * size_seaweed + 1, 4 * size_seaweed + 2 });
		shape.connectivity.push_back({ 4 * size_seaweed + 3, 4 * size_seaweed, 4 * size_seaweed + 2 });

		normal(shape.position, shape.connectivity, shape.normal);
		seaweeds.push_back(shape);
	}

	return seaweeds;
}

mesh_drawable scene_model::create_glass(vec3 p0, vec3 u0, vec3 u1, vec3 u2)
{

	mesh_drawable glass = mesh_drawable(mesh_primitive_parallelepiped(p0, u0, u1, u2));
	glass.uniform.color_alpha = 0.05f;
	glass.uniform.color = vec3(0.5f, 0.5f, 0.5f);
	glass.uniform.shading.specular = 0.8f;
	glass.uniform.shading.ambiant = 0.8f;
	glass.uniform.shading.diffuse = 0.8f;
	return glass;
}

void scene_model::setup_data(std::map<std::string, GLuint>& shaders, scene_structure& scene, gui_structure& gui)
{

	// Initial position of the camera
	scene.camera.translation = { 0.0f, -0.5f, 0.0f };
	scene.camera.scale = 5.0f;
	gui.show_frame_camera = false;

	//Water

	int n_water = 201;
	water = create_water(n_water);

	water_draw = water;
	water_draw.shader = shaders["mesh"];
	water_draw.uniform.color_alpha = 0.3f;
	water_draw.uniform.shading.specular = 0.6f;
	water_draw.uniform.shading.ambiant = 0.7f;
	water_draw.uniform.shading.diffuse = 0.5f;

	//Sand

	int n_sand = 201;
	sand = create_sand(n_sand);

	sand_draw = sand;
	sand_draw.shader = shaders["mesh"];
	sand_draw.uniform.color_alpha = 1.0f;
	sand_draw.uniform.shading.specular = 0.0f;
	sand_draw.uniform.shading.ambiant = 0.5f;
	sand_draw.uniform.shading.diffuse = 0.5f;
	sand_draw.uniform.transform.scaling = 0.999f;

	//Glass
	glass_back = create_glass({ -1.05, -0.05, -1.05 }, { 2.1, 0, 0 }, { 0, 1.25, 0 }, { 0, 0, 0.05 });
	glass_front = create_glass({ -1.05, -0.05, 1 }, { 2.1, 0, 0 }, { 0, 1.25, 0 }, { 0, 0, 0.05 });
	glass_left = create_glass({ -1.05, -0.05, -1.05 }, { 0.05, 0, 0 }, { 0, 1.25, 0 }, { 0, 0, 2.1 });
	glass_right = create_glass({ 1, -0.05, -1.05 }, { 0.05, 0, 0 }, { 0, 1.25, 0 }, { 0, 0, 2.1 });
	glass_bottom = create_glass({ -1.05, -0.05, -1.05 }, { 2.1, 0, 0 }, { 0, 0.05, 0 }, { 0, 0, 2.1 });


	//Seaweeds

	seaweeds = create_seaweeds(nb_seaweeds, size_seaweed);

	// Init seaweeds data (speed, force)

	initial_positions_seaweeds.resize(nb_seaweeds);
	for (int i = 0; i < nb_seaweeds; i++)
	{
		for (int j = 0; j < 4 * (size_seaweed + 1); j++)
		{
			initial_positions_seaweeds[i].push_back(seaweeds[i].position[j]);
		}
	}

	speeds_seaweeds.resize(nb_seaweeds);
	for (int i = 0; i < nb_seaweeds; i++)
	{
		speeds_seaweeds[i].resize(4 * (size_seaweed + 1));
		speeds_seaweeds[i].fill({ 0, 0, 0 });
	}

	forces_seaweeds.resize(nb_seaweeds);
	for (int i = 0; i < nb_seaweeds; i++)
	{
		forces_seaweeds[i].resize(4 * (size_seaweed + 1));
		forces_seaweeds[i].fill({ 0, 0, 0 });
	}

	for (int i = 0; i < nb_seaweeds; i++)
	{
		seaweeds_draw.push_back(seaweeds[i]);
		seaweeds_draw[i].uniform.color = vec3(0.2f, 0.7f, 0.2f);
		seaweeds_draw[i].uniform.color_alpha = 1.0f;
		seaweeds_draw[i].uniform.shading.specular = 0.2f;
		seaweeds_draw[i].uniform.shading.ambiant = 0.4f;
		seaweeds_draw[i].uniform.shading.diffuse = 0.2f;
	}

	//Textures

	texture_sand = create_texture_gpu(image_load_png("scenes/sources/default/animation/assets/rocks.png"));

	sphere = mesh_drawable(mesh_primitive_sphere(1.0f));
	sphere.shader = shaders["mesh"];

	//Fishes

	//initialize position and speed of the fishes
	init_fishes();

	// create the hierarchical structures of the fishes
	for (int i = 0; i < N_fishes; i++) {
		hierarchy_fishes.push_back(create_fish(shaders));
	}
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

static vec3 spring_force(const vec3& pi, const vec3& pj, float L0, float K)
{
	const vec3 p = pi - pj;
	float L = norm(p);
	const vec3 u = normalize(p);

	const vec3 F = -K * (L - L0) * u;
	return F;
}

void scene_model::compute_forces_seaweeds(int number_of_seaweeds, int size_seaweed)
{

	// Get simuation parameters
	const float K = 3.5f; // /!\ Divergence if too high ! No more than 4.f
	const float Krigid = 20.f;
	const float m = 0.01f / (4 * size_seaweed);
	const float wind_force = 5.5f;
	const float rand_wind_prop = 0.3;
	const float mu = 0.0f;

	// Reset

	for (size_t k = 0; k < number_of_seaweeds; k++)
	{
		for (size_t i = 1; i < 4 * (size_seaweed + 1); i++)
		{
			forces_seaweeds[k][i] = vec3(0, 0, 0);
		}
	}

	// Drag

	for (size_t k = 0; k < number_of_seaweeds; k++)
	{
		for (size_t i = 1; i < 4 * (size_seaweed + 1); i++)
		{
			forces_seaweeds[k][i] -= mu * (speeds_seaweeds[k][i]);
		}
	}

	// Wind

	std::default_random_engine generator;
	std::uniform_real_distribution<double> distribution(-1, 1);

	for (size_t k = 0; k < number_of_seaweeds; k++)
	{

		float L0 = (initial_positions_seaweeds[k][4 * size_seaweed] - initial_positions_seaweeds[k][0])[1] / size_seaweed;
		float L1 = (initial_positions_seaweeds[k][2] - initial_positions_seaweeds[k][0])[0];

		vec3 rand_wind = vec3(cos(10 * k) * rand_wind_prop * wind_force * distribution(generator), 0, sin(10 * k) * rand_wind_prop * wind_force * distribution(generator));

		for (size_t i = 1; i < 4 * (size_seaweed + 1); i++)
		{
			float rand_phase = 0.4 * distribution(generator);
			vec3 wind = vec3(cos(10 * k + timer.t + rand_phase) * wind_force, 0.3 * wind_force, sin(10 * k + timer.t) * wind_force + rand_phase) + rand_wind;
			vec3 wind_u = normalize(wind);
			vec3& n = seaweeds[k].normal[i];
			float wind_magnitude = std::abs(dot(wind, n));
			vec3 f = wind_magnitude * wind_u;

			forces_seaweeds[k][i] += f;
		}
	}

	// Springs

	for (size_t k = 0; k < number_of_seaweeds; k++)
	{
		float L0 = (initial_positions_seaweeds[k][4 * size_seaweed] - initial_positions_seaweeds[k][0])[1] / size_seaweed;
		float L1 = (initial_positions_seaweeds[k][2] - initial_positions_seaweeds[k][0])[0];
		float L2 = sqrt(L0 * L0 + L1 * L1);
		float L3 = sqrt(2) * L1;

		for (size_t i = 1; i < size_seaweed; i++)
		{
			//Direct Neighbour springs
			forces_seaweeds[k][4 * i] += spring_force(seaweeds[k].position[4 * i], seaweeds[k].position[4 * (i - 1)], L0, K / L0);
			forces_seaweeds[k][4 * i] += spring_force(seaweeds[k].position[4 * i], seaweeds[k].position[4 * (i + 1)], L0, K / L0);
			forces_seaweeds[k][4 * i] += spring_force(seaweeds[k].position[4 * i], seaweeds[k].position[4 * i + 1], L1, K / L1);
			forces_seaweeds[k][4 * i] += spring_force(seaweeds[k].position[4 * i], seaweeds[k].position[4 * i + 3], L1, K / L1);

			forces_seaweeds[k][4 * i + 1] += spring_force(seaweeds[k].position[4 * i + 1], seaweeds[k].position[4 * (i - 1) + 1], L0, K / L0);
			forces_seaweeds[k][4 * i + 1] += spring_force(seaweeds[k].position[4 * i + 1], seaweeds[k].position[4 * (i + 1) + 1], L0, K / L0);
			forces_seaweeds[k][4 * i + 1] += spring_force(seaweeds[k].position[4 * i + 1], seaweeds[k].position[4 * i], L1, K / L1);
			forces_seaweeds[k][4 * i + 1] += spring_force(seaweeds[k].position[4 * i + 1], seaweeds[k].position[4 * i + 2], L1, K / L1);

			forces_seaweeds[k][4 * i + 2] += spring_force(seaweeds[k].position[4 * i + 2], seaweeds[k].position[4 * (i - 1) + 2], L0, K / L0);
			forces_seaweeds[k][4 * i + 2] += spring_force(seaweeds[k].position[4 * i + 2], seaweeds[k].position[4 * (i + 1) + 2], L0, K / L0);
			forces_seaweeds[k][4 * i + 2] += spring_force(seaweeds[k].position[4 * i + 2], seaweeds[k].position[4 * i + 1], L1, K / L1);
			forces_seaweeds[k][4 * i + 2] += spring_force(seaweeds[k].position[4 * i + 2], seaweeds[k].position[4 * i + 3], L1, K / L1);

			forces_seaweeds[k][4 * i + 3] += spring_force(seaweeds[k].position[4 * i + 3], seaweeds[k].position[4 * (i - 1) + 3], L0, K / L0);
			forces_seaweeds[k][4 * i + 3] += spring_force(seaweeds[k].position[4 * i + 3], seaweeds[k].position[4 * (i + 1) + 3], L0, K / L0);
			forces_seaweeds[k][4 * i + 3] += spring_force(seaweeds[k].position[4 * i + 3], seaweeds[k].position[4 * i + 2], L1, K / L1);
			forces_seaweeds[k][4 * i + 3] += spring_force(seaweeds[k].position[4 * i + 3], seaweeds[k].position[4 * i], L1, K / L1);

			//Diagonal Springs (vertically speaking)

			forces_seaweeds[k][4 * i] += spring_force(seaweeds[k].position[4 * i], seaweeds[k].position[4 * (i - 1) + 1], L2, K / L2);
			forces_seaweeds[k][4 * i] += spring_force(seaweeds[k].position[4 * i], seaweeds[k].position[4 * (i - 1) + 3], L2, K / L2);
			forces_seaweeds[k][4 * i] += spring_force(seaweeds[k].position[4 * i], seaweeds[k].position[4 * (i + 1) + 1], L2, K / L2);
			forces_seaweeds[k][4 * i] += spring_force(seaweeds[k].position[4 * i], seaweeds[k].position[4 * (i + 1) + 3], L2, K / L2);

			forces_seaweeds[k][4 * i + 1] += spring_force(seaweeds[k].position[4 * i + 1], seaweeds[k].position[4 * (i - 1)], L2, K / L2);
			forces_seaweeds[k][4 * i + 1] += spring_force(seaweeds[k].position[4 * i + 1], seaweeds[k].position[4 * (i - 1) + 2], L2, K / L2);
			forces_seaweeds[k][4 * i + 1] += spring_force(seaweeds[k].position[4 * i + 1], seaweeds[k].position[4 * (i + 1)], L2, K / L2);
			forces_seaweeds[k][4 * i + 1] += spring_force(seaweeds[k].position[4 * i + 1], seaweeds[k].position[4 * (i + 1) + 2], L2, K / L2);

			forces_seaweeds[k][4 * i + 2] += spring_force(seaweeds[k].position[4 * i + 2], seaweeds[k].position[4 * (i - 1) + 1], L2, K / L2);
			forces_seaweeds[k][4 * i + 2] += spring_force(seaweeds[k].position[4 * i + 2], seaweeds[k].position[4 * (i - 1) + 3], L2, K / L2);
			forces_seaweeds[k][4 * i + 2] += spring_force(seaweeds[k].position[4 * i + 2], seaweeds[k].position[4 * (i + 1) + 1], L2, K / L2);
			forces_seaweeds[k][4 * i + 2] += spring_force(seaweeds[k].position[4 * i + 2], seaweeds[k].position[4 * (i + 1) + 3], L2, K / L2);

			forces_seaweeds[k][4 * i + 3] += spring_force(seaweeds[k].position[4 * i + 3], seaweeds[k].position[4 * (i - 1) + 2], L2, K / L2);
			forces_seaweeds[k][4 * i + 3] += spring_force(seaweeds[k].position[4 * i + 3], seaweeds[k].position[4 * (i - 1)], L2, K / L2);
			forces_seaweeds[k][4 * i + 3] += spring_force(seaweeds[k].position[4 * i + 3], seaweeds[k].position[4 * (i + 1) + 2], L2, K / L2);
			forces_seaweeds[k][4 * i + 3] += spring_force(seaweeds[k].position[4 * i + 3], seaweeds[k].position[4 * (i + 1)], L2, K / L2);

			//Diagonal Springs (horizontally speaking)

			forces_seaweeds[k][4 * i] += spring_force(seaweeds[k].position[4 * i], seaweeds[k].position[4 * i + 2], L3, K / L3);
			forces_seaweeds[k][4 * i + 1] += spring_force(seaweeds[k].position[4 * i + 1], seaweeds[k].position[4 * i + 3], L3, K / L3);
			forces_seaweeds[k][4 * i + 2] += spring_force(seaweeds[k].position[4 * i + 2], seaweeds[k].position[4 * i], L3, K / L3);
			forces_seaweeds[k][4 * i + 3] += spring_force(seaweeds[k].position[4 * i + 3], seaweeds[k].position[4 * i + 1], L3, K / L3);

			//Rigidity springs (to force the seaweed to stay at initial position)

			forces_seaweeds[k][4 * i] += spring_force(seaweeds[k].position[4 * i], initial_positions_seaweeds[k][4 * i], 0, Krigid / i);
			forces_seaweeds[k][4 * i + 1] += spring_force(seaweeds[k].position[4 * i + 1], initial_positions_seaweeds[k][4 * i + 1], 0, Krigid / i);
			forces_seaweeds[k][4 * i + 2] += spring_force(seaweeds[k].position[4 * i + 2], initial_positions_seaweeds[k][4 * i + 2], 0, Krigid / i);
			forces_seaweeds[k][4 * i + 3] += spring_force(seaweeds[k].position[4 * i + 3], initial_positions_seaweeds[k][4 * i + 3], 0, Krigid / i);
		}

		for (size_t i = 1; i < size_seaweed; i++)
		{
			//2-Distance springs

			if (i > 1)
			{
				forces_seaweeds[k][4 * i] += spring_force(seaweeds[k].position[4 * i], seaweeds[k].position[4 * (i - 2)], 2 * L0, K / (L0 * 2));
				forces_seaweeds[k][4 * i + 1] += spring_force(seaweeds[k].position[4 * i + 1], seaweeds[k].position[4 * (i - 2) + 1], 2 * L0, K / (L0 * 2));
				forces_seaweeds[k][4 * i + 2] += spring_force(seaweeds[k].position[4 * i + 2], seaweeds[k].position[4 * (i - 2) + 2], 2 * L0, K / (L0 * 2));
				forces_seaweeds[k][4 * i + 3] += spring_force(seaweeds[k].position[4 * i + 3], seaweeds[k].position[4 * (i - 2) + 3], 2 * L0, K / (L0 * 2));
			}
			if (i < size_seaweed - 1)
			{
				forces_seaweeds[k][4 * i] += spring_force(seaweeds[k].position[4 * i], seaweeds[k].position[4 * (i + 2)], 2 * L0, K / (L0 * 2));
				forces_seaweeds[k][4 * i + 1] += spring_force(seaweeds[k].position[4 * i + 1], seaweeds[k].position[4 * (i + 2) + 1], 2 * L0, K / (L0 * 2));
				forces_seaweeds[k][4 * i + 2] += spring_force(seaweeds[k].position[4 * i + 2], seaweeds[k].position[4 * (i + 2) + 2], 2 * L0, K / (L0 * 2));
				forces_seaweeds[k][4 * i + 3] += spring_force(seaweeds[k].position[4 * i + 3], seaweeds[k].position[4 * (i + 2) + 3], 2 * L0, K / (L0 * 2));
			}
		}

		//Last layer

		forces_seaweeds[k][4 * size_seaweed] += spring_force(seaweeds[k].position[4 * size_seaweed], seaweeds[k].position[4 * (size_seaweed - 1)], L0, K / L0);
		forces_seaweeds[k][4 * size_seaweed] += spring_force(seaweeds[k].position[4 * size_seaweed], seaweeds[k].position[4 * size_seaweed + 1], L1, K / L1);
		forces_seaweeds[k][4 * size_seaweed] += spring_force(seaweeds[k].position[4 * size_seaweed], seaweeds[k].position[4 * size_seaweed + 3], L1, K / L1);

		forces_seaweeds[k][4 * size_seaweed + 1] += spring_force(seaweeds[k].position[4 * size_seaweed + 1], seaweeds[k].position[4 * (size_seaweed - 1) + 1], L0, K / L0);
		forces_seaweeds[k][4 * size_seaweed + 1] += spring_force(seaweeds[k].position[4 * size_seaweed + 1], seaweeds[k].position[4 * size_seaweed], L1, K / L1);
		forces_seaweeds[k][4 * size_seaweed + 1] += spring_force(seaweeds[k].position[4 * size_seaweed + 1], seaweeds[k].position[4 * size_seaweed + 2], L1, K / L1);

		forces_seaweeds[k][4 * size_seaweed + 2] += spring_force(seaweeds[k].position[4 * size_seaweed + 2], seaweeds[k].position[4 * (size_seaweed - 1) + 2], L0, K / L0);
		forces_seaweeds[k][4 * size_seaweed + 2] += spring_force(seaweeds[k].position[4 * size_seaweed + 2], seaweeds[k].position[4 * size_seaweed + 1], L1, K / L1);
		forces_seaweeds[k][4 * size_seaweed + 2] += spring_force(seaweeds[k].position[4 * size_seaweed + 2], seaweeds[k].position[4 * size_seaweed + 3], L1, K / L1);

		forces_seaweeds[k][4 * size_seaweed + 3] += spring_force(seaweeds[k].position[4 * size_seaweed + 3], seaweeds[k].position[4 * (size_seaweed - 1) + 3], L0, K / L0);
		forces_seaweeds[k][4 * size_seaweed + 3] += spring_force(seaweeds[k].position[4 * size_seaweed + 3], seaweeds[k].position[4 * size_seaweed + 2], L1, K / L1);
		forces_seaweeds[k][4 * size_seaweed + 3] += spring_force(seaweeds[k].position[4 * size_seaweed + 3], seaweeds[k].position[4 * size_seaweed], L1, K / L1);

		//Diagonal Springs (vertically speaking)

		forces_seaweeds[k][4 * size_seaweed] += spring_force(seaweeds[k].position[4 * size_seaweed], seaweeds[k].position[4 * (size_seaweed - 1) + 1], L2, K / L2);
		forces_seaweeds[k][4 * size_seaweed] += spring_force(seaweeds[k].position[4 * size_seaweed], seaweeds[k].position[4 * (size_seaweed - 1) + 3], L2, K / L2);

		forces_seaweeds[k][4 * size_seaweed + 1] += spring_force(seaweeds[k].position[4 * size_seaweed + 1], seaweeds[k].position[4 * (size_seaweed - 1)], L2, K / L2);
		forces_seaweeds[k][4 * size_seaweed + 1] += spring_force(seaweeds[k].position[4 * size_seaweed + 1], seaweeds[k].position[4 * (size_seaweed - 1) + 2], L2, K / L2);

		forces_seaweeds[k][4 * size_seaweed + 2] += spring_force(seaweeds[k].position[4 * size_seaweed + 2], seaweeds[k].position[4 * (size_seaweed - 1) + 1], L2, K / L2);
		forces_seaweeds[k][4 * size_seaweed + 2] += spring_force(seaweeds[k].position[4 * size_seaweed + 2], seaweeds[k].position[4 * (size_seaweed - 1) + 3], L2, K / L2);

		forces_seaweeds[k][4 * size_seaweed + 3] += spring_force(seaweeds[k].position[4 * size_seaweed + 3], seaweeds[k].position[4 * (size_seaweed - 1) + 2], L2, K / L2);
		forces_seaweeds[k][4 * size_seaweed + 3] += spring_force(seaweeds[k].position[4 * size_seaweed + 3], seaweeds[k].position[4 * (size_seaweed - 1)], L2, K / L2);

		//Diagonal Springs (horizontally speaking)

		forces_seaweeds[k][4 * size_seaweed] += spring_force(seaweeds[k].position[4 * size_seaweed], seaweeds[k].position[4 * size_seaweed + 2], L3, K / L3);
		forces_seaweeds[k][4 * size_seaweed + 1] += spring_force(seaweeds[k].position[4 * size_seaweed + 1], seaweeds[k].position[4 * size_seaweed + 3], L3, K / L3);
		forces_seaweeds[k][4 * size_seaweed + 2] += spring_force(seaweeds[k].position[4 * size_seaweed + 2], seaweeds[k].position[4 * size_seaweed], L3, K / L3);
		forces_seaweeds[k][4 * size_seaweed + 3] += spring_force(seaweeds[k].position[4 * size_seaweed + 3], seaweeds[k].position[4 * size_seaweed + 1], L3, K / L3);

		//Rigidity springs (to force the seaweed to stay at initial position)

		forces_seaweeds[k][4 * size_seaweed] += spring_force(seaweeds[k].position[4 * size_seaweed], initial_positions_seaweeds[k][4 * size_seaweed], 0, Krigid / size_seaweed);
		forces_seaweeds[k][4 * size_seaweed + 1] += spring_force(seaweeds[k].position[4 * size_seaweed + 1], initial_positions_seaweeds[k][4 * size_seaweed + 1], 0, Krigid / size_seaweed);
		forces_seaweeds[k][4 * size_seaweed + 2] += spring_force(seaweeds[k].position[4 * size_seaweed + 2], initial_positions_seaweeds[k][4 * size_seaweed + 2], 0, Krigid / size_seaweed);
		forces_seaweeds[k][4 * size_seaweed + 3] += spring_force(seaweeds[k].position[4 * size_seaweed + 3], initial_positions_seaweeds[k][4 * size_seaweed + 3], 0, Krigid / size_seaweed);
	}
}

void scene_model::update_seaweeds(int nb_seaweeds, int size_seaweed, float dt)
{

	const float m = 0.1f / (4 * size_seaweed);

	for (size_t i = 0; i < nb_seaweeds; ++i)
	{

		for (size_t j = 12; j < 4 * (size_seaweed + 1); ++j)
		{
			vec3& p = seaweeds[i].position[j];
			vec3& v = speeds_seaweeds[i][j];
			const vec3& f = forces_seaweeds[i][j];

			v = v + dt * f / m;
			if (norm(v) > 0.01)
			{
				normalize(v);
				v *= 0.01;
			}
			p = p + dt * v;
		}

		normal(seaweeds[i].position, seaweeds[i].connectivity, seaweeds[i].normal, true);

		seaweeds_draw[i].update_position(seaweeds[i].position);
		seaweeds_draw[i].update_normal(seaweeds[i].normal);
	}
}

void scene_model::frame_draw(std::map<std::string, GLuint>& shaders, scene_structure& scene, gui_structure&)
{

	float dt = 0.005f;
	if (!(timer.update() > 0))
		dt = 0;

	if (f)
	{

		display_fish(shaders,scene);
	}

	compute_forces_seaweeds(nb_seaweeds, size_seaweed);
	update_seaweeds(nb_seaweeds, size_seaweed, dt);

	if (sw)
	{
		for (int i = 0; i < nb_seaweeds; i++)
		{

			draw(seaweeds_draw[i], scene.camera, shaders["mesh"]);
			if (gui_scene.wireframe)
				draw(seaweeds_draw[i], scene.camera, shaders["wireframe"]);
		}
	}

	if (s)
	{

		glBindTexture(GL_TEXTURE_2D, texture_sand);

		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_MIRRORED_REPEAT);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_MIRRORED_REPEAT);

		draw(sand_draw, scene.camera);
		if (gui_scene.wireframe)
			draw(sand_draw, scene.camera, shaders["wireframe"]);
	}

	glBindTexture(GL_TEXTURE_2D, scene.texture_white);

	set_gui();
	update_water(dt);
	compute_time_step_fishes(dt);
	timer.t += dt;
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
	ImGui::Checkbox("Fish", &f);
	ImGui::Checkbox("Sand", &s);
	ImGui::Checkbox("Seaweeds", &sw);
}


void scene_model::init_fishes() {
	for (size_t k = 0; k < N_fishes; k++) {
		particle_structure fish;


		// Initial random position
		float x = rand_interval(0, 0.5);
		float y = rand_interval(0, 0.5);
		float z = rand_interval(0, 0.5);
		fish.p = vec3(x, y, z);
		std::cout << fish.p << std::endl;

		// Initial random speed
		float vx = rand_interval(0, 1);
		float vy = rand_interval(0, 1);
		float vz = rand_interval(0, 1);
		fish.v = vec3(vx, vy, vz);

		// no forces at the beginning
		fish.f = vec3(0, 0, 0);

		fishes.push_back(fish);
	}


}

//function f for the first Boids model we tried (not used in the final version)
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
	return alpha1 * exp(-k * (pow((x - x0) / x0, 2))) - alpha2 * exp(-x / x0);

}
void scene_model::compute_time_step_fishes(float dt)
{
	// Set forces
	const size_t N = fishes.size();

	//first tried model (not used)
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

	//definition of the zone for each rule of the Boids model
	float rsep = 0.5;
	float ralign = 1.0;
	float rcohez = 2.0;

	//maximal speed
	float vmax = 5.0f;

	//add boid speed for each fish
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
		fish.f = -alpha_frot * fish.v;

		//attraction to the center
		vec3 center = { 0,0.5,0 };
		float raideur = 0.1;
		float l0 = 0.1f;
		fish.f += -raideur * (norm(p - center) - l0) * (fish.p - center);

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

		//coefficients for the importance of each effect of the Boids model
		float alpha_sep = 0.05;
		float alpha_align = 0.01;
		float alpha_cohez = 0.01;


		//add to the speed
		v += alpha_sep * vsep + alpha_align * valign + alpha_align * vcohez;

		//set a maximal norm to avoid divergence
		if (norm(v) > vmax) {
			v = v * vmax / norm(v);
		}

		//random distrubance
		float alpha_alea = 0.05;
		float vx = rand_interval(-1, 1);
		float vy = rand_interval(-1, 1);
		float vz = rand_interval(-1, 1);
		fish.v += alpha_alea * vec3(vx, vy, vz);

		//repuslive force on the edges of the aquarium
		float drepuls = 0.2f;
		float frepuls = 100.0f;
		vec3 force_repuls = { 0,0,0 };
		if (p.x > 1 - drepuls) {
			force_repuls += vec3(-1, 0, 0);
		}
		if (p.x < -1 + drepuls) {
			force_repuls += vec3(1, 0, 0);
		}
		if (p.y > 1 - drepuls) {
			force_repuls += vec3(0, -1, 0);
		}
		if (p.y < drepuls) {
			force_repuls += vec3(0, 1, 0);
		}
		if (p.z > 1 - drepuls) {
			force_repuls += vec3(0, 0, -1);
		}
		if (p.z < -1 + drepuls) {
			force_repuls += vec3(0, 0, 1);
		}

		//zero if the fish moves away from the edge
		float valeur = -dot(force_repuls, fish.v) / norm(fish.v);

		if (valeur > 0) {
			fish.f += force_repuls * valeur * frepuls;
		}
	}

	// Integrate position and speed of particles through time
	for (size_t k = 0; k < N; ++k) {
		particle_structure& fish = fishes[k];
		vec3& v = fish.v;
		vec3& p = fish.p;
		vec3 const& f = fish.f;

		v = v + dt * f;
		p = p + dt * v;
	}

	//Borders :
	float atte = 0.0f;
	float r = 0.035f;
	float h_sand = 0.1f;
	for (size_t k = 0; k < N; ++k) {
		particle_structure& fish = fishes[k];
		vec3& p = fish.p;
		vec3& v = fish.v;
		if (p.x > 1 - r) {
			p.x = 1 - r;
			v.x = -atte * v.x;
		}
		if (p.x < -1 + r) {
			p.x = -1 + r;
			v.x = -atte * v.x;
		}
		if (p.y > 1 - r) {
			p.y = 1 - r;
			v.y = -atte * v.y;
		}
		if (p.y < r+h_sand) {
			p.y = r+h_sand;
			v.y = -atte * v.y;
		}
		if (p.z > 1 - r) {
			p.z = 1 - r;
			v.z = -atte * v.z;
		}
		if (p.z < -1 + r) {
			p.z = -1 + r;
			v.z = -atte * v.z;
		}

	}


}


//create the structure hierarchy of the fish
hierarchy_mesh_drawable scene_model::create_fish(std::map<std::string, GLuint>& shaders)
{
	hierarchy_mesh_drawable hierarchy;
	float r = rand_interval(0.02f,0.03f); //to have fishes of different sizes
	float vert = rand_interval(0.0f, 1.0f); //to have fishes of different colors from yellow to red
	vec3 c = { 1.0f,vert,0 };

	mesh fish_body = mesh_primitive_sphere(r);
	for (size_t k = 0; k < fish_body.position.size(); ++k)
	{
		fish_body.position[k].x *= 1.5f;
		fish_body.position[k].z *= 0.75f;
	}
	fish_body.fill_color_uniform(c);
	mesh_drawable body = fish_body;
	hierarchy.add(body, "body", "root", { 0, 0, 0 });


	mesh fish_tail = mesh_primitive_cone(2 * r / 3, { -2 * r, 0, 0 }, { -r, 0.0f, -0.0f });
	fish_tail.fill_color_uniform(c);
	hierarchy.add(fish_tail, "tail", "body", { 0.0f, 0.0f, 0.0f });

	mesh fish_eyeR = mesh_primitive_sphere(r / 5);
	fish_eyeR.fill_color_uniform({ 0,0,0 });
	hierarchy.add(fish_eyeR, "fish_eyeR", "body", { r,r / 2,r / 2 });

	mesh fish_eyeL = mesh_primitive_sphere(r / 5);
	fish_eyeL.fill_color_uniform({ 0,0,0 });
	mesh_drawable fish_eyeL_drawable = fish_eyeL;
	hierarchy.add(fish_eyeL_drawable, "fish_eyeL", "body", { r, r / 2,-r / 2 });

	mesh finL = mesh_primitive_quad({ -2 * r / 3, 0.0f, 0 }, { -2 * r / 3, 0.0f, -r / 4 }, { r / 2, 0.0f, -3 * r / 4 }, { r / 2,0.0f,0 });
	finL.fill_color_uniform(c);
	hierarchy.add(finL, "finL", "body", { 0.0f, 0.0f, -r / 2 });


	mesh finR = mesh_primitive_quad({ -2 * r / 3, 0.0f,0 }, { -2 * r / 3, 0.0f, r / 4 }, { r / 2, 0.0f, 3 * r / 4 }, { r / 2,0.0f,0 });
	finR.fill_color_uniform(c);
	hierarchy.add(finR, "finR", "body", { 0.0f,0.0f,r / 2 });
	hierarchy.update_local_to_global_coordinates();
	hierarchy.set_shader_for_all_elements(shaders["mesh"]);
	return hierarchy;
}

//to animate and display the fishes
void scene_model::display_fish(std::map<std::string, GLuint>& shaders, scene_structure& scene)
{


	const size_t N = fishes.size();
	for (size_t k = 0; k < N; ++k)
	{
		const particle_structure& part = fishes[k];

		const float theta = std::cos(3 * 3.14f * timer.t + k); //k : phase to have a different movement between the fishes

		//rotation of the lateral fin
		affine_transform& local_finL = hierarchy_fishes[k]["finL"].transform;
		local_finL.rotation = rotation_from_axis_angle_mat3({ 1, 0, 0 }, theta);

		affine_transform& local_finR = hierarchy_fishes[k]["finR"].transform;
		local_finR.rotation = rotation_from_axis_angle_mat3({ 1, 0, 0 }, -theta);

		//rotation of the tail
		affine_transform& local_tail = hierarchy_fishes[k]["tail"].transform;
		local_tail.rotation = rotation_from_axis_angle_mat3({ 0, 1, 0 }, theta / 3.0f) * rotation_from_axis_angle_mat3({ 1, 0, 0 }, theta);

		//translation of the body
		affine_transform& local_body = hierarchy_fishes[k]["body"].transform;
		local_body.translation = part.p;

		//global orientation of the body in the direction of the speed
		if (norm(part.v) > 0) {
			vec3 vit = part.v / norm(part.v);
			local_body.rotation = rotation_between_vector_mat3({ 1,0,0 }, vit);
		}
		hierarchy_fishes[k].update_local_to_global_coordinates();
		draw(hierarchy_fishes[k], scene.camera);
	}


}

#endif