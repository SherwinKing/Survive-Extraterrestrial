#include "WalkMesh.hpp"

#include "read_write_chunk.hpp"

#include <glm/gtx/norm.hpp>
#include <glm/gtx/string_cast.hpp>

#include <iostream>
#include <fstream>
#include <algorithm>
#include <string>

WalkMesh::WalkMesh(std::vector< glm::vec3 > const &vertices_, std::vector< glm::vec3 > const &normals_, std::vector< glm::uvec3 > const &triangles_)
	: vertices(vertices_), normals(normals_), triangles(triangles_) {

	//construct next_vertex map (maps each edge to the next vertex in the triangle):
	next_vertex.reserve(triangles.size()*3);
	auto do_next = [this](uint32_t a, uint32_t b, uint32_t c) {
		auto ret = next_vertex.insert(std::make_pair(glm::uvec2(a,b), c));
		assert(ret.second);
	};
	for (auto const &tri : triangles) {
		do_next(tri.x, tri.y, tri.z);
		do_next(tri.y, tri.z, tri.x);
		do_next(tri.z, tri.x, tri.y);
	}

	// //DEBUG: are vertex normals consistent with geometric normals?
	// for (auto const &tri : triangles) {
	// 	glm::vec3 const &a = vertices[tri.x];
	// 	glm::vec3 const &b = vertices[tri.y];
	// 	glm::vec3 const &c = vertices[tri.z];
	// 	glm::vec3 out = glm::normalize(glm::cross(b-a, c-a));

	// 	float da = glm::dot(out, normals[tri.x]);
	// 	float db = glm::dot(out, normals[tri.y]);
	// 	float dc = glm::dot(out, normals[tri.z]);

	// 	assert(da > 0.1f && db > 0.1f && dc > 0.1f);
	// }
}

//project pt to the plane of triangle a,b,c and return the barycentric weights of the projected point:
glm::vec3 barycentric_weights(glm::vec3 const &a, glm::vec3 const &b, glm::vec3 const &c, glm::vec3 const &pt) {
	//TODO: implement!
	//compute barycentric weights here!
	//ref: https://gamedev.stackexchange.com/questions/23743/whats-the-most-efficient-way-to-find-barycentric-coordinates
	glm::vec3 ab = b - a;
	glm::vec3 ac = c - a;
	glm::vec3 ap = pt - a;
	
	float d_ab_ab = glm::dot(ab, ab);
	float d_ab_ac = glm::dot(ab, ac);
	float d_ac_ac = glm::dot(ac, ac);
	float d_ap_ab = glm::dot(ap, ab);
	float d_ap_ac = glm::dot(ap, ac);
	float denom = d_ab_ab * d_ac_ac - d_ab_ac * d_ab_ac;
	float inv_denom = 1 / denom;


	float v = (d_ac_ac * d_ap_ab - d_ab_ac * d_ap_ac) * inv_denom;
	float w = (d_ab_ab * d_ap_ac - d_ab_ac * d_ap_ab) * inv_denom;
	float u = 1.0f - v - w;
	return glm::vec3(u, v, w);
	
	// my previous algo.
	/*
	{
	glm::vec3 ab = b - a;
	glm::vec3 ac = c - a;
	glm::vec3 bc = c - b;
	
	// Norm of the cross product is the triangle area.
	glm::vec3 cross_ab_ac = glm::cross(ab, ac);
	
	
	
	glm::vec3 ap = pt - a;
	glm::vec3 bp = pt - b;
	
	
	
	float bary_a = glm::cross(bc, bp)[0] / cross_ab_ac[0];
	float bary_b = glm::cross(ap, ac)[0] / cross_ab_ac[0];
	float bary_c = glm::cross(ab, bp)[0] / cross_ab_ac[0];

	return glm::vec3(bary_a, bary_b, bary_c);
	}
	*/
}

WalkPoint WalkMesh::nearest_walk_point(glm::vec3 const &world_point) const {
	assert(!triangles.empty() && "Cannot start on an empty walkmesh");

	WalkPoint closest;
	float closest_dis2 = std::numeric_limits< float >::infinity();

	for (auto const &tri : triangles) {
		//find closest point on triangle:

		glm::vec3 const &a = vertices[tri.x];
		glm::vec3 const &b = vertices[tri.y];
		glm::vec3 const &c = vertices[tri.z];

		//get barycentric coordinates of closest point in the plane of (a,b,c):
		glm::vec3 coords = barycentric_weights(a,b,c, world_point);

		//is that point inside the triangle?
		if (coords.x >= 0.0f && coords.y >= 0.0f && coords.z >= 0.0f) {
			//yes, point is inside triangle.
			float dis2 = glm::length2(world_point - to_world_point(WalkPoint(tri, coords)));
			if (dis2 < closest_dis2) {
				closest_dis2 = dis2;
				closest.indices = tri;
				closest.weights = coords;
			}
		} else {
			//check triangle vertices and edges:
			auto check_edge = [&world_point, &closest, &closest_dis2, this](uint32_t ai, uint32_t bi, uint32_t ci) {
				glm::vec3 const &a = vertices[ai];
				glm::vec3 const &b = vertices[bi];

				//find closest point on line segment ab:
				float along = glm::dot(world_point-a, b-a);
				float max = glm::dot(b-a, b-a);
				glm::vec3 pt;
				glm::vec3 coords;
				if (along < 0.0f) {
					pt = a;
					coords = glm::vec3(1.0f, 0.0f, 0.0f);
				} else if (along > max) {
					pt = b;
					coords = glm::vec3(0.0f, 1.0f, 0.0f);
				} else {
					float amt = along / max;
					pt = glm::mix(a, b, amt);
					coords = glm::vec3(1.0f - amt, amt, 0.0f);
				}

				float dis2 = glm::length2(world_point - pt);
				if (dis2 < closest_dis2) {
					closest_dis2 = dis2;
					closest.indices = glm::uvec3(ai, bi, ci);
					closest.weights = coords;
				}
			};
			check_edge(tri.x, tri.y, tri.z);
			check_edge(tri.y, tri.z, tri.x);
			check_edge(tri.z, tri.x, tri.y);
		}
	}
	assert(closest.indices.x < vertices.size());
	assert(closest.indices.y < vertices.size());
	assert(closest.indices.z < vertices.size());
	return closest;
}


void WalkMesh::walk_in_triangle(WalkPoint const &start, glm::vec3 const &step, WalkPoint *end_, float *time_) const {
	assert(end_);
	auto &end = *end_;

	assert(time_);
	auto &time = *time_;

	glm::vec3 const &a = vertices[start.indices.x];
	glm::vec3 const &b = vertices[start.indices.y];
	glm::vec3 const &c = vertices[start.indices.z];
	
	// imagining taking the step from a
	glm::vec3 a_step = a + step;
	// glm::vec3 is barycentric coord for a.
	glm::vec3 barycentric_velocity = barycentric_weights(a, b, c, a_step) - glm::vec3(1, 0, 0);
	//assert(barycentric_velocity[0] + barycentric_velocity[1] + barycentric_velocity[2] == 0);

	//TODO: check when/if this velocity pushes start.weights into an edge
	glm::vec3 possible_end = start.weights + barycentric_velocity;
	float back_step_size = 0;
	float current_back_step_size;
	for (int i = 0; i < 3; i++) {
		if (possible_end[i] < 0) {
			current_back_step_size = (possible_end[i] / barycentric_velocity[i]);
			glm::vec3 back_bary_step = -barycentric_velocity * current_back_step_size;
			possible_end += back_bary_step;
			possible_end[i] = 0;
			back_step_size += current_back_step_size;
		}
	}
	end.indices = start.indices;
	end.weights = possible_end;
	time = 1 - back_step_size;


	//Remember: our convention is that when a WalkPoint is on an edge,
	// then wp.weights.z == 0.0f (so will likely need to re-order the indices)
	
	// Check if follow convention. If not, re-order.
	if (back_step_size > 0) {
		if (end.weights.z != 0.0f) {
			if (end.weights.x == 0.0f) {
				end.indices = glm::uvec3(end.indices.y, end.indices.z, end.indices.x);
				end.weights = glm::vec3(end.weights.y, end.weights.z, end.weights.x);
			} else if (end.weights.y == 0.0f) {
				end.indices = glm::uvec3(end.indices.z, end.indices.x, end.indices.y);
				end.weights = glm::vec3(end.weights.z, end.weights.x, end.weights.y);
			}
		}
	}
}

bool WalkMesh::cross_edge(WalkPoint const &start, WalkPoint *end_, glm::quat *rotation_) const {
	assert(start.weights.z == 0.0f);
	assert(start.indices.x <= vertices.size() && start.indices.y <= vertices.size() && start.indices.z <= vertices.size());
	assert(end_);
	auto &end = *end_;
	assert(rotation_);
	auto &rotation = *rotation_;
//!todo{

	//TODO: check if edge (start.indices.x, start.indices.y) has a triangle on the other side:
	//  hint: remember 'next_vertex'!
	auto twin_tri_vertex_it = next_vertex.find(glm::uvec2(start.indices.y, start.indices.x));
	if (twin_tri_vertex_it == next_vertex.end())
		return false;
	
	//TODO: if there is another triangle:
	//  TODO: set end's weights and indicies on that triangle:
	end.indices = {start.indices.y, start.indices.x, twin_tri_vertex_it->second};
	end.weights = {start.weights.y, start.weights.x, 0};
	

	//  TODO: compute rotation that takes starting triangle's normal to ending triangle's normal:
	//  hint: look up 'glm::rotation' in the glm/gtx/quaternion.hpp header
	// ref: https://stackoverflow.com/questions/1171849/finding-quaternion-representing-the-rotation-from-one-vector-to-another
	glm::vec3 s_xy = vertices[start.indices.y] - vertices[start.indices.x];
	glm::vec3 s_xz = vertices[start.indices.z] - vertices[start.indices.x];
	glm::vec3 e_xy = vertices[end.indices.y] - vertices[end.indices.x];
	glm::vec3 e_xz = vertices[end.indices.z] - vertices[end.indices.x];
	
	glm::vec3 norm_s_unnormalized = glm::cross(s_xy, s_xz);
	glm::vec3 norm_e_unnormalized = glm::cross(e_xy, e_xz);
	glm::vec3 norm_cross_unnormalized = glm::cross(norm_s_unnormalized, norm_e_unnormalized);
	float w = glm::length(norm_s_unnormalized) * glm::length(norm_e_unnormalized) + glm::dot(norm_s_unnormalized, norm_e_unnormalized);

	rotation = glm::normalize(glm::quat(w, norm_cross_unnormalized[0], norm_cross_unnormalized[1], norm_cross_unnormalized[2])); //identity quat (wxyz init order)

	//return 'true' if there was another triangle, 'false' otherwise:
	return true;
}


WalkMeshes::WalkMeshes(std::string const &filename) {
	std::ifstream file(filename, std::ios::binary);

	std::vector< glm::vec3 > vertices;
	read_chunk(file, "p...", &vertices);

	std::vector< glm::vec3 > normals;
	read_chunk(file, "n...", &normals);

	std::vector< glm::uvec3 > triangles;
	read_chunk(file, "tri0", &triangles);

	std::vector< char > names;
	read_chunk(file, "str0", &names);

	struct IndexEntry {
		uint32_t name_begin, name_end;
		uint32_t vertex_begin, vertex_end;
		uint32_t triangle_begin, triangle_end;
	};

	std::vector< IndexEntry > index;
	read_chunk(file, "idxA", &index);

	if (file.peek() != EOF) {
		std::cerr << "WARNING: trailing data in walkmesh file '" << filename << "'" << std::endl;
	}

	//-----------------

	if (vertices.size() != normals.size()) {
		throw std::runtime_error("Mis-matched position and normal sizes in '" + filename + "'");
	}

	for (auto const &e : index) {
		if (!(e.name_begin <= e.name_end && e.name_end <= names.size())) {
			throw std::runtime_error("Invalid name indices in index of '" + filename + "'");
		}
		if (!(e.vertex_begin <= e.vertex_end && e.vertex_end <= vertices.size())) {
			throw std::runtime_error("Invalid vertex indices in index of '" + filename + "'");
		}
		if (!(e.triangle_begin <= e.triangle_end && e.triangle_end <= triangles.size())) {
			throw std::runtime_error("Invalid triangle indices in index of '" + filename + "'");
		}

		//copy vertices/normals:
		std::vector< glm::vec3 > wm_vertices(vertices.begin() + e.vertex_begin, vertices.begin() + e.vertex_end);
		std::vector< glm::vec3 > wm_normals(normals.begin() + e.vertex_begin, normals.begin() + e.vertex_end);

		//remap triangles:
		std::vector< glm::uvec3 > wm_triangles; wm_triangles.reserve(e.triangle_end - e.triangle_begin);
		for (uint32_t ti = e.triangle_begin; ti != e.triangle_end; ++ti) {
			if (!( (e.vertex_begin <= triangles[ti].x && triangles[ti].x < e.vertex_end)
			    && (e.vertex_begin <= triangles[ti].y && triangles[ti].y < e.vertex_end)
			    && (e.vertex_begin <= triangles[ti].z && triangles[ti].z < e.vertex_end) )) {
				throw std::runtime_error("Invalid triangle in '" + filename + "'");
			}
			wm_triangles.emplace_back(
				triangles[ti].x - e.vertex_begin,
				triangles[ti].y - e.vertex_begin,
				triangles[ti].z - e.vertex_begin
			);
		}
		
		std::string name(names.begin() + e.name_begin, names.begin() + e.name_end);

		auto ret = meshes.emplace(name, WalkMesh(wm_vertices, wm_normals, wm_triangles));
		if (!ret.second) {
			throw std::runtime_error("WalkMesh with duplicated name '" + name + "' in '" + filename + "'");
		}

	}
}

WalkMesh const &WalkMeshes::lookup(std::string const &name) const {
	auto f = meshes.find(name);
	if (f == meshes.end()) {
		throw std::runtime_error("WalkMesh with name '" + name + "' not found.");
	}
	return f->second;
}
