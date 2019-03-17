#include <fstream>
#include <sstream>
#include <map>
#include <iostream>
#include <algorithm>
#include <cmath>

#include "WingEdge.h"


bool WingEdge::loadOBJfile(std::string filename)
{
	std::vector<vec3> vertices;
	std::vector<vec3> normals;
	std::vector<std::vector<int>> faces;

	std::ifstream ifstr(filename);
	if (!ifstr)
		return false;

	vertices.clear();
	faces.clear();
	normals.clear();

	// load vertices and faces
	std::string line;
	while (std::getline(ifstr, line))
	{
		if (line.length() < 1 || line[0] == '#')
			continue;

		else if (line[0] == 'v')
		{
			// load a vertex
			line.erase(line.begin());
			std::istringstream iss(line);
			vec3 v;
			iss >> v.x >> v.y >> v.z;
			vertices.push_back(v);

		}
		else if (line[0] == 'f')
		{
			// load a face
			line.erase(line.begin());
			faces.resize(faces.size() + 1);
			std::istringstream iss(line);
			int idx;
			while (iss >> idx)
				faces.back().push_back(idx - 1);
		}
	}

	// Estimate normals of vertices with the weights defined by angles around them
	normals = std::vector<vec3>(vertices.size(), vec3(0, 0, 0));
	for (int i = 0; i < faces.size(); i++)
	{
		int vnum = faces[i].size();
		for (int j = 0; j < vnum; j++) {
			vec3 v = vertices[faces[i][j]];
			vec3 v_pre = vertices[faces[i][(j - 1 + vnum) % vnum]];
			vec3 v_nex = vertices[faces[i][(j + 1) % vnum]];
			double weight = (v_pre - v).vectorDegreeBetween(v_nex - v);
			normals[faces[i][j]] = normals[faces[i][j]] + ((v - v_pre) ^ (v_nex - v)).normalize() * weight;
		}
	}
	for (int i = 0; i < normals.size(); ++i)
		normals[i].normalize();

	normalizeEdge(vertices, normals, faces);

	// convert the mesh to winged-edge form
	convertOBJToWingedEdgeMesh(vertices, normals, faces);

	return true;
}

void WingEdge::normalizeEdge(std::vector<vec3> &vertices, std::vector<vec3> &normals, std::vector<std::vector<int>> &faces)
{
	// find bbox of the mesh
	vec3 min(0, 0, 0);
	vec3 max(0, 0, 0);
	for (int i = 0; i < vertices.size(); i++)
	{
		vec3 v = vertices[i];
		if (v.x < min.x) min.x = v.x;
		if (v.y < min.y) min.y = v.y;
		if (v.z < min.z) min.z = v.z;

		if (v.x > max.x) max.x = v.x;
		if (v.y > max.y) max.y = v.y;
		if (v.z > max.z) max.z = v.z;
	}

	center = (max + min) / 2;
	vec3 dimlengths = max - min;

	maxdimlength = dimlengths.x;
	if (dimlengths.y > maxdimlength) maxdimlength = dimlengths.y;
	if (dimlengths.z > maxdimlength) maxdimlength = dimlengths.z;

	// normalize bounding box of the vertices
	for (int i = 0; i < vertices.size(); i++)
		vertices[i] = (vertices[i] - center) / maxdimlength;
}

void WingEdge::convertOBJToWingedEdgeMesh(std::vector<vec3> &vertices, std::vector<vec3> &normals, std::vector<std::vector<int>> &faces)
{
	// create wvertex and save to vertexlist
	vertexList.clear();
	for (int i = 0; i < vertices.size(); i++)
		vertexList.push_back(new Wvertex(vertices[i], normals[i], i));

	// create wedge and save to edge
	edgeList.clear();
	for (int i = 0; i < faces.size() * 3; i++)
		edgeList.push_back(new Wedge());

	// create wface and save to face
	faceList.clear();
	for (int i = 0; i < faces.size(); i++)
		faceList.push_back(new Wface());

	int vnum = faces[0].size();

	// use edgetable to record indices of start and end vertex, edge index
	int m = vertices.size();
	std::vector<std::vector<int>> edgeTable(m, std::vector<int>(m, -1));

	for (int i = 0; i < faces.size(); i++)
	{
		Wface *f = faceList[i];

		for (int v = 0; v < vnum; v++)
		{
			int v0 = faces[i][v];
			Wvertex *v_end = vertexList[v0];

			int v1 = faces[i][(v + 1) % vnum];
			Wvertex *v_start = vertexList[v1];

			Wedge *edge = edgeList[vnum*i + v];
			edge->v_start = v_start;
			edge->v_end = v_end;

			Wedge *edge_nex = edgeList[vnum*i + (v - 1 + vnum) % vnum];
			edge->e_right_nex = edge_nex;

			Wedge *edge_pre = edgeList[vnum*i + (v + 1) % vnum];
			edge->e_right_pre = edge_pre;

			edgeTable[v1][v0] = vnum * i + v;
			edge->f_right = f;

			f->oneEdge = edge;
			v_start->oneEdge = edge;
			// v_end->oneEdge = edge; 
		}
	}

	// save left according the table
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < m; j++)
		{
			int idx = edgeTable[i][j];
			if (idx != -1)
			{
				Wedge *edge = edgeList[idx];
				int flipIdx = edgeTable[j][i];
				Wedge *edgeFlip = edgeList[flipIdx];

				edge->e_left_pre = edgeFlip->e_right_pre;
				edge->e_left_nex = edgeFlip->e_right_nex;

				edge->f_left = edgeFlip->f_right;
			}
		}
	}

}

std::vector<std::vector<Wvertex *>> WingEdge::extractVerticesOfFaces()
{
	// the vertices index
	std::vector<std::vector<Wvertex *>> faces(faceList.size());

	for (int i = 0; i < faceList.size(); i++)
	{
		Wface *f = faceList[i];
		Wedge *e0 = f->oneEdge;
		Wedge *edge = e0;

		// find all the edges related to the face
		do
		{
			Wface *ft = edge->f_right;

			if (edge->f_right == f)
				edge = edge->e_right_pre;

			Wvertex *vt = edge->v_start;
			faces[i].push_back(edge->v_start);
		} while (edge != e0);
	}

	////// update norms ////////////
	for (int i = 0; i < faces.size(); i++)
	{
		int vnum = faces[i].size();
		for (int j = 0; j < vnum; j++) {
			vec3 v = faces[i][j]->position;
			vec3 v_pre = faces[i][(j - 1 + vnum) % vnum]->position;
			vec3 v_nex = faces[i][(j + 1) % vnum]->position;
			double weight = (v_pre - v).vectorDegreeBetween(v_nex - v);
			faces[i][j]->norm = faces[i][j]->norm + ((v - v_pre) ^ (v_nex - v)).normalize() * weight;
		}
	}

	for (int i = 0; i < faces.size(); i++)
		for (int j = 0; j < faces[i].size(); j++)
			faces[i][j]->norm.normalize();

	return faces;
}

bool WingEdge::saveToOBJfile(std::string fname)
{
	// check whether we could create the file
	std::ofstream ofstr(fname);
	if (!ofstr)
		return false;

	// extract the face list from winged-edge mesh
	std::vector<std::vector<Wvertex *>> faces = extractVerticesOfFaces();

	// export vertices
	for (int i = 0; i < vertexList.size(); i++)
	{
		vec3 v = (vertexList[i]->position* maxdimlength) + center;
		ofstr << "v " << v.x
			<< " " << v.y
			<< " " << v.z << std::endl;
	}

	int idx;

	// export faces
	for (int i = 0; i < faces.size(); i++)
	{
		ofstr << "f";

		for (int j = 0; j < faces[i].size(); j++)
		{
			idx = faces[i][j]->idx;
			ofstr << " " << idx + 1;
		}
		ofstr << std::endl;
	}

	return true;
}

// for map
bool operator<(const int2 a, const int2 b) {
	if (a.x < b.x) return true;
	if (a.x > b.x) return false;
	if (a.y < b.y) return true;
	return false;
}
bool operator==(const int2 a, const int2 b) {
	return (a.x == b.x) && (a.y == b.y);
}


bool WingEdge::decimation(int k, int target)
{
	bool result = false;
	computeQuadrics();

	while (target--)
	{
		result = multipleChoice(k);
	}
	
	return result;
}

void WingEdge::computeQuadrics()
{
	for (int i = 0; i < faceList.size(); i++)
	{
		Wface *f = faceList[i];
		Wedge *e = f->oneEdge;
		vec3 v1 = e->v_start->position;
		vec3 v2 = e->v_end->position;
		vec3 v3 = e->e_right_nex->v_end->position;
		vec3 cross = ((v3 - v1) ^ (v2 - v1)).normalize();

		float a = cross.x;
		float b = cross.y;
		float c = cross.z;
		float d = -(v1.x*a + v1.y*b + v1.z*c);

		// for each vertex of the face
		Wedge *edge = e;
		int j = 0;
		do {
			Wvertex *v = e->v_start;
			Eigen::Matrix4f k;
			k << a*a, a*b, a*c, a*d,
				a*b, b*b, b*c, b*d,
				a*c, b*c, c*c, c*d,
				a*d, b*d, c*d, d*d;
			v->Q = v->Q + k;
			e = e->e_right_nex;
		} while (edge != e);

	}
}

bool WingEdge::multipleChoice(int k)
{
	// select k random candidates
	std::vector<Wedge *> collapseCandidates;

	int i = k;
	int j = 0;
	while (i > 0)
	{
		Wedge *e = edgeList[(rand() % edgeList.size())];
		j++;
		if (checkCandidate(e))
		{
			collapseCandidates.push_back(e);
			i--;
		}
		if (j > edgeList.size())
			break;
	}

	//find the edge with the least error
	float currentError = INFINITY;
	Wedge *currentEdge = NULL;
	vec3 currentPosition;
	Eigen::Vector4f position;

	for (int i = 0; i < collapseCandidates.size(); i++)
	{
		Wedge *edge = collapseCandidates[i];

		float error = computeError(edge, position);
		//std::cout << error << std::endl;
		if (error < currentError)
		{
			currentError = error;
			currentEdge = edge;
			currentPosition.x = position(0);
			currentPosition.y = position(1);
			currentPosition.z = position(2);
		}
	}
	
	// collapse the edge
	if (currentEdge != NULL)
	{
		collapseEdge(currentEdge, currentPosition);
		return true;
	}
	else
	{
		std::cout << "can not find an edge to collapse" << std::endl;
		return false;
	}
}

bool WingEdge::checkCandidate(Wedge *edge)
{
	Wvertex *va = edge->e_right_nex->v_end;
	Wvertex *vb = edge->e_left_nex->v_end;

	if (getValence(va) == 3 || getValence(vb) == 3)
		return false;
	return true;
}

int WingEdge::getValence(Wvertex *vertex)
{
	int v = 0;

	Wedge *e = vertex->oneEdge;
	Wedge *e0 = e;

	do
	{
		e = e->e_left_nex;
		v += 1;
	} while (e != e0);

	return v;
}

float WingEdge::computeError(Wedge *edge, Eigen::Vector4f &position)
{
	Eigen::Matrix4f Q = edge->v_start->Q + edge->v_end->Q;
	Eigen::Vector4f v;
	Eigen::Matrix4f temp = Q;
	temp.row(3) = Eigen::Vector4f(0.0, 0.0, 0.0, 1.0); //
	Eigen::Matrix4f inverseQ = temp.inverse();
	float error;

	// check if Q is inversible
	if (std::isnan(inverseQ(0, 0)))
	{
		// use midpoint or endpoints
		Eigen::Vector4f v0 = Eigen::Vector4f(edge->v_start->position.x, edge->v_start->position.y, edge->v_start->position.z, 1.0);
		Eigen::Vector4f v1 = Eigen::Vector4f(edge->v_end->position.x, edge->v_end->position.y, edge->v_end->position.z, 1.0);
		Eigen::Vector4f vm = Eigen::Vector4f((v0(0) + v1(0))*0.5, (v0(1) + v1(1))*0.5, (v0(2) + v1(2))*0.5, 1.0);

		float e0 = v0.transpose() * Q * v0;
		float e1 = v1.transpose() * Q * v1;
		float em = vm.transpose() * Q * vm;

		if (e0 < e1 && e0 < em)
		{
			error = e0;
			position = v0;
		}
		else if (e1 < e0 && e1 < em)
		{
			error = e1;
			position = v1;
		}
		else
		{
			error = em;
			position = vm;
		}
	}
	else
	{
		v = inverseQ * Eigen::Vector4f(0.0, 0.0, 0.0, 1.0);
		error = v.transpose() * Q * v;
		position = v;
	}

	return error;
}

void WingEdge::collapseEdge(Wedge *edge, vec3 position)
{
	// "new" vertex
	Wvertex *vnew = edge->v_end;
	vnew->position = position;
	vnew->Q += edge->v_start->Q;
	vnew->oneEdge = edge->e_left_nex->e_left_nex;

	//// update vertex
	Wedge *e = edge->e_left_nex->e_left_nex;

	int count = getValence(edge->v_start);

	if (count >= 3)
	{
		count = 0;
		// update vertex idx
		while (e != edge)
		{
			e->v_start = vnew;
			e->e_right_pre->v_end = vnew;
			e = e->e_left_nex;
			count += 1;
		}
	}

	//// update four edge relations
	Wedge *edge_up_right = edge->e_right_pre->e_left_pre->e_right_nex; // flip(rightp)

	Wedge *edge_up_left = edge->e_right_nex->e_left_pre->e_right_nex; // flip(rightn)

	Wedge *edge_down_right = edge->e_left_nex->e_left_pre->e_right_nex; // flip(left next)

	Wedge *edge_down_left = edge->e_left_pre->e_left_pre->e_right_nex; // flip(left pre)

	// up right
	edge_up_right->e_left_nex = edge_up_left->e_right_nex;
	edge_up_right->e_left_pre = edge_up_left->e_right_pre;

	// up left
	edge_up_left->e_left_nex = edge_up_right->e_right_nex;
	edge_up_left->e_left_pre = edge_up_right->e_right_pre;

	// down right
	edge_down_right->e_left_nex = edge_down_left->e_right_nex;
	edge_down_right->e_left_pre = edge_down_left->e_right_pre;

	// down left
	edge_down_left->e_left_nex = edge_down_right->e_right_nex;
	edge_down_left->e_left_pre = edge_down_right->e_right_pre;

	// update vertex->one edge
	edge_up_left->v_start->oneEdge = edge_up_left;
	edge_down_right->v_start->oneEdge = edge_down_right;

	//// update edge face
	edge_up_right->f_left = edge_up_left->f_right;
	edge_up_left->f_left = edge_up_right->f_right;

	edge_down_right->f_left = edge_down_left->f_right;
	edge_down_left->f_left = edge_down_right->f_right;

	//// detele vertex, edge, face and erase from the list
	Wvertex *vs = edge->v_start;
	Wface *fup = edge->f_right;
	Wface *fdown = edge->f_left;

	// delete edges
	edge->e_right_nex->v_start = vs;
	edge->e_left_pre->v_start = vs;

	std::vector<Wedge *>::iterator eit;
	for (eit = edgeList.begin(); eit != edgeList.end();)
	{
		if ((*eit)->v_start == vs || (*eit)->v_end == vs)
		{
			delete *eit;
			eit = edgeList.erase(eit);
		}
		else
			eit++;
	}

	// delete faces
	std::vector<Wface *>::iterator fit;
	for (fit = faceList.begin(); fit != faceList.end();)
	{
		if ((fup && (*fit) == fup) || (fdown && (*fit) == fdown))
		{
			delete *fit;
			fit = faceList.erase(fit);
		}
		else
			fit++;
	}

	// delete vertex
	int idxs = vs->idx;
	for (int k = idxs+1; k < vertexList.size(); k++)
		vertexList[k]->idx -= 1;

	delete vertexList[idxs];
	vertexList.erase(vertexList.begin() + idxs);
}