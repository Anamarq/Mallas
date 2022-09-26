/*
 * meshBoundary.cpp
 *
 * Written by Jose Miguel Espadero <josemiguel.espadero@urjc.es>
 *
 * This code is written as material for the FMF class of the
 * Master Universitario en Informatica Grafica, Juegos y Realidad Virtual.
 * Its purpose is to be didactic and easy to understand, not hard optimized.
 *
 * This file computes the list of internal and external edges of a mesh
 * and writes the boundary (list of external edges). Then create a copy
 * of the mesh over a colorMesh and assign the red color to the vertex
 * in the boundary.
 * 
 * //TODO: Fill-in your name and email
 * Name of alumn: Ana Márquez Moncada
 * Email of alumn: anamarq98@gmail.com
 * Year: 2021
 * 
 */

#ifdef _MSC_VER
#pragma warning(error: 4101)
#endif

#define _CRT_NONSTDC_NO_DEPRECATE
#include <iostream>
#include <cmath>
#include <SimpleMesh.hpp>
#include <ColorMesh.hpp>
#include <chrono>
#include <queue>
#include<set>

#include<algorithm>
#include<vector>
using namespace std::chrono;

/// Update the contents of externalEdges and internalEdges 
void updateEdgeLists(const SimpleMesh &mesh, 
                     std::vector<SimpleEdge> &externalEdges,
                     std::vector<SimpleEdge> &internalEdges )
{

    //TODO 2.1: Implement the body of the updateEdgeLists() method
	std::set<SimpleEdge>externalAux;
	std::set<SimpleEdge>internalAux;
	for (int t = 0; t < mesh.numTriangles(); ++t) {
		for (int e = 0; e < 3; ++e) {
			//Si la arista pertenece a externalEdges o internalEdges
			auto ed = mesh.triangles[t].edges()[e];
			//Searches the container for an element equivalent to val and returns an iterator to it if found, otherwise it returns an iterator to set::end.
			auto indexExt=externalAux.find(ed);
			auto indexInt= internalAux.find(ed);
			if ((indexExt != externalAux.end()) || (indexInt != internalAux.end())) {
				cout << "Arista repetida, la malla no es mainfold" << endl;
			}
			//Si la arista inversa está en externalEdges
			auto indexInverse = externalAux.find(ed.reversed());
			if ( indexInverse!= externalAux.end()){
				//Borrar arista inversa de externalEdges
				externalAux.erase(indexInverse );
				//Insertar la arista en internalEdges
				internalAux.insert(ed.reversed());
			}		
			else {
				//Insertar arista en externalEdges
				externalAux.insert(ed);
			}
		}
	}
	externalEdges.insert(externalEdges.end(), externalAux.begin(), externalAux.end());
	internalEdges.insert(internalEdges.end(), internalAux.begin(), internalAux.end());

	if (externalEdges.size() > 0) {
		//El vértice “a” de la primera arista será el de menor coordenada X de todos los vértices de
		//frontera.En caso de empate, aquel de ellos que tenga menor coordenada Y.En caso de
		//	empate, aquel de ellos que tenga menor coordenada Z
		float menorX = mesh.coordinates[externalEdges[0].a].X;
		float menorY = mesh.coordinates[externalEdges[0].a].Y;
		float menorZ = mesh.coordinates[externalEdges[0].a].Z;
		int idMenor = 0;
		for (int i = 1; i < externalEdges.size(); ++i) {
			float actualX = mesh.coordinates[externalEdges[i].a].X;
			float actualY = mesh.coordinates[externalEdges[i].a].Y;
			float actualZ = mesh.coordinates[externalEdges[i].a].Z;
			if (actualX < menorX) {
				menorX = actualX;
				idMenor = i;
				menorY = actualY;
			}
			else if (actualX == menorX) {
				if (actualY < menorY) {
					menorX = actualX;
					menorY = actualY;
					menorZ = actualZ;
					idMenor = i;
				}
				else if (actualY == menorY) {
					if (actualZ < menorZ) {
						menorX = actualX;
						menorY = actualY;
						menorZ = actualZ;
						idMenor = i;
					}
				}
			}
		}
		//intercabio posición (primer elemento)
		SimpleEdge aux = externalEdges[idMenor];
		externalEdges[idMenor] = externalEdges[0];
		externalEdges[0] = aux;
		//El resto de aristas forman una (o varias) listas circulares, es decir, el índice “a” de cada
		//arista corresponderá con el índice “b” de la arista anterior(ver ejemplo más adelante)
		for (int i = 1; i < externalEdges.size(); i++) {
			// Buscar la posicion donde esta el vertice, donde coinciden el vertice b en pos i-1 con el vertice a en la pos k
			int k = i;
			while ((externalEdges[i - 1].b != externalEdges[k].a) && (k < externalEdges.size()))
			{
				++k;
			}

			if (k == externalEdges.size())  //Esto pasa cuando es una malla de múltiples fronteras, ya que no habrá continuidad en todas 
				k = i;                     //las aristas externas, por lo que se deja en la misma posición
			// intercambiar pos
			SimpleEdge aux = externalEdges[i];
			externalEdges[i] = externalEdges[k];
			externalEdges[k] = aux;
		}
	}
    //throw ("updateEdgeLists has to be implemented as exercise");
    //END TODO 2.1

}//void updateEdgeLists()



//Algoritmo de Dijkstra, devuelve el coste minimo de los vertices a los externos
void Dijkstra(SimpleMesh &mesh, const int v_ini, const int v_fin,std::vector<double>&distancia, std::vector<std::vector<int>>&neighbour)
{
	priority_queue<pair<double, int>,vector<pair<double,int>>,greater<pair<double,int>>>pq;  // coste y nodo
	vector<int> padre(mesh.numVertex());      
	
	pq.push(make_pair(0,v_ini));  //Añadir nodo inicial a la cola (coste 0)
	//mientras que no este la cola vacía...
	while (!pq.empty()) {
		double c = pq.top().first;  //coste del nodo u
		int u = pq.top().second; //nodo con mínimo coste
		pq.pop(); //eliminar de la cola 

		//actualiza con el camino más corto, mirando los vecinos directos del vértice u
		for (int v = 0; v <neighbour[u].size(); v++) {
			int vv = neighbour[u][v];  //vecino del vértice u
			//Si hay un camino más corto para vv a través de u
			if (distancia[vv] > distancia[u] + mesh.distance(u, vv)) {
				distancia[vv] = distancia[u] + mesh.distance(u, vv); //Nueva distancia
				padre[vv] = u;
				pq.push(make_pair( distancia[vv],vv ) ); //insertar vv en pq
			}
		}
	}
}


int main (int argc, char *argv[])
{
    try
    {
        // Set default input mesh filename
        //std::string filename("mallas/16Triangles.off");
       // std::string filename("mallas/mannequin2.ply");
       // std::string filename("mallas/knot-hole.ply");
        //std::string filename("c:/Nefertiti.990kv.ply");
		std::string filename("mallas/mask2.ply");
		//std::string filename("mallas/bunny.ply");
		//std::string filename("c:/angel_kneeling.150kv.ply");
        if (argc > 1)
            filename = std::string(argv[1]);

        ///////////////////////////////////////////////////////////////////////
        //Read a mesh and write a mesh
        SimpleMesh mesh;        
        cout << "Loading file " << filename << endl;
        mesh.readFile(filename, false);

        cout << "Num vertex: " << mesh.numVertex() << " Num triangles: " << mesh.numTriangles()
             << " Unreferenced vertex: " << mesh.checkUnreferencedVertex() << endl;

        //Init time measures used for profiling
        high_resolution_clock::time_point clock0 = high_resolution_clock::now();

        ///////////////////////////////////////////////////////////////////////

        //Compute the list of external and internal edges
        //Implement the body of the updateEdgeLists() method .
        std::vector<SimpleEdge> externalEdges;
        std::vector<SimpleEdge> internalEdges;
        updateEdgeLists(mesh, externalEdges, internalEdges);

 
        cout << "Done updateEdgeLists() " << duration<float>(high_resolution_clock::now() - clock0).count() << " seconds" << endl;
        cout << externalEdges.size() << " boundary edges" << endl;
        cout << internalEdges.size() << " internal edges" << endl;

        //Write external boundary and compute boundary length
        double boundaryLength = 0.0;
        cout << "Edges in the boundary:" << endl;
        for (const SimpleEdge &e : externalEdges)
        {
            cout << "[" << e.a << "->" << e.b << "] ";
            boundaryLength += mesh.edgeLength(e);
        }
        cout << endl;
        cout << "Boundary length: " << boundaryLength << endl;

        //Compute the Euler Characteristic for this mesh
        //For a conected mesh, the value means:
        //2 -> The mesh is a closed surface with no holes (sphere-like topology)
        //1 -> The mesh has one hole (disk-like topology)
        //0 -> The has one handle (torus-like topology)
        //-2 -> The mesh has two holes (tube-like topology)
        int eulerCharacteristic = mesh.numVertex() + mesh.numTriangles()
                                - externalEdges.size() - internalEdges.size();
        cout << "Euler Characteristic: " << eulerCharacteristic << endl;


        //Vector to store the distance from a vertex to the nearest vertex in the boundary
        std::vector<double> boundDist;
        //Index of the vertex with max distance to boundary
        unsigned deepestVertex = 0xDEADC0DE;

        //TODO 2.2:
        //Compute the distance from each vertex to the nearest vertex in the boundary (euclidean distance measured along edges)
        //and store it in the boundDist vector. 
        //Store in deepestVertex the index of the vertex with the maximum distance to boundary
		boundDist.resize(mesh.numVertex()); //vector de distancias
		std::vector<int> externalVertices; //vector que almacena los vértices externos
		//externalVertices.reserve(externalEdges.size());
		for (int i = 0; i < externalEdges.size(); ++i) {
			externalVertices.push_back(externalEdges[i].a);  //añado el vertice a de cada arista externa
		}
		//inicializar vector de distancia 
		for (int i = 0; i < mesh.numVertex(); i++) {
			boundDist[i] = INFINITY;
		}
		//las distancias desde vertices externos es 0 (ya es frontera)
		for (int i = 0; i < externalVertices.size(); i++) {
			boundDist[externalVertices[i]] = 0;
		}
		std::vector<std::vector<int>> neighbour;
		mesh.computeNeighbours(neighbour);//calcula los vecinos de cada vértice


		double costeMin;
		double costeMax = 0;  //coste del vértice con mayor coste (entre todos los vértices)
		//bucle que recorre todos los vértices externos
		for (int i = 0; i < externalVertices.size(); ++i) {
			costeMin = INFINITY;  //guarda el coste minimo de un vertice a uno externo (un solo vértice)
			Dijkstra(mesh, externalVertices[i], -1, boundDist,neighbour);  //Caminos más cortos de los vertices a los externos				
		}

		for (int ii = 0; ii< boundDist.size(); ++ii) {
			if (boundDist[ii] > costeMax) {
				costeMax = boundDist[ii];
				deepestVertex = ii;
			}
		}
						
        //END TODO 2.2
        cout << "Done boundDist() " << duration<float>(high_resolution_clock::now() - clock0).count() << " seconds" << endl;

        //Dump the index and distance for the deepestVertex
        cout << "maxDistance to boundary: " << boundDist[deepestVertex] << " at vertex: " << deepestVertex << endl;

        //Dump distances to boundary (for small meshes) for debugging
        if (boundDist.size() < 40)
        {
            cout << "Distances to boundary: " << endl;
            for (size_t i=0; i < boundDist.size(); i++)
                 cout << "vertex: " << i << " : " << boundDist[i] << endl;
            cout << endl;
        }

        /* Dump path from deepestVertex to nearest boundary vertex
        cout << "Shortest path from " << deepestVertex<< " to boundary:" << endl;
        int i = deepestVertex;
        while (i != -1)
        {
            cout << i << " -> ";
            i = parent[i];
        }
        cout << endl;
        //Dump path */

		
			


        std::string outputFilename="output_boundary.ply";
        cout << "Saving output to " << outputFilename << endl;

        //TODO 2.3:
        //Create a color mesh where the color of each vertex shows its distance to boundary
        //Save the color mesh to a PLY file named "output_boundary.ply"
        //See meshColor.cpp or meshColor2.cpp to see an how-to example

		//Prepare a copy of the input mesh on a ColorMesh
		ColorMesh outputMesh;
		outputMesh.coordinates = mesh.coordinates;
		outputMesh.triangles = mesh.triangles;
		outputMesh.colors.resize(mesh.numVertex());

		//Iterate over the vertex of the mesh and set the color
		//of each vertex using its distance to basePoint.
		if (externalEdges.size() > 0) {             //Si esto pasa es que no hay agujero
			for (size_t i = 0; i < mesh.numVertex(); i++)
			{
				outputMesh.colors[i].setTemperature(float(costeMax - boundDist[i]), 0, float(costeMax));
			}
		}
		
		///////////////////////////////////////////////////////////////////////
		//Save result to a file in .ply format

		outputMesh.writeFilePLY(outputFilename);



	

        //END TODO 2.3
        cout << "Done colorize by distance to boundary " << duration<float>(high_resolution_clock::now() - clock0).count() << " seconds" << endl;

        //Visualize the file with an external viewer
#ifdef WIN32
        //string viewcmd = "\"C:\\Program Files (x86)\\VCG\\MeshLab\\meshlab.exe\"";
        string viewcmd = "C:/meshlab/meshlab.exe";
#else
        string viewcmd = "meshlab >/dev/null 2>&1 ";
#endif
        string cmd = viewcmd+" "+outputFilename;
        cout << "Executing external command: " << cmd << endl;
        return system(cmd.c_str());
    }

    catch (const string &str) { std::cerr << "EXCEPTION: " << str << std::endl; }
    catch (const char *str) { std::cerr << "EXCEPTION: " << str << std::endl; }
    catch (std::exception& e)    { std::cerr << "EXCEPTION: " << e.what() << std::endl;  }
    catch (...) { std::cerr << "EXCEPTION (unknow)" << std::endl; }

#ifdef WIN32
    cout << "Press Return to end the program" <<endl;
    cin.get();
#else
#endif

    return 0;
}

