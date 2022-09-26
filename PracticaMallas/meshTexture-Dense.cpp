/*
 * meshTexture-Dense.cpp
 *
 * Written by Jose Miguel Espadero <josemiguel.espadero@urjc.es>
 *
 * This code is written as material for the FMF class of the
 * Master Universitario en Informatica Grafica, Juegos y Realidad Virtual.
 * Its purpose is to be didactic and easy to understand, not hard optimized.
 * 
 * //TODO: Fill-in your name and email
 * Name of alumn: Ana Márquez Moncada
 * Email of alumn: anamarq98@gmail.com
 * Year: 2021
 * 
 */

//This file is the simplest solution to the meshTexture exercise.
//Use a dense matrix and does not reduce the problem, so take long time to execute
//for files as mannequin2.ply

#define _CRT_NONSTDC_NO_DEPRECATE
#include <iostream>
#include <iomanip>
#include <chrono>
#include<set>
using namespace std::chrono;

#include <SimpleMesh.hpp>
#include <TextureMesh.hpp>

//Check if Eigen  (a standar Matrix Library) is included. If not,
//you can get it at http://eigen.tuxfamily.org

// <Eigen/Dense> is the module for dense (traditional) matrix and vector.
// You can get a quick reference for using Eigen dense objects at:
// http://eigen.tuxfamily.org/dox/group__QuickRefPage.html
#include <Eigen/Core>
#include <Eigen/LU>

/// Update the contents of externalEdges and internalEdges 
void updateEdgeLists(const SimpleMesh &mesh, 
                     std::vector<SimpleEdge> &externalEdges,
                     std::vector<SimpleEdge> &internalEdges )
{
    //TODO 2.1: Copy the body of the updateEdgeLists() method from meshBoundary


	std::set<SimpleEdge>externalAux;
	std::set<SimpleEdge>internalAux;
	for (int t = 0; t < mesh.numTriangles(); ++t) {
		for (int e = 0; e < 3; ++e) {
			//Si la arista pertenece a externalEdges o internalEdges
			auto ed = mesh.triangles[t].edges()[e];
			//Searches the container for an element equivalent to val and returns an iterator to it if found, otherwise it returns an iterator to set::end.
			auto indexExt = externalAux.find(ed);
			auto indexInt = internalAux.find(ed);
			if ((indexExt != externalAux.end()) || (indexInt != internalAux.end())) {
				cout << "Arista repetida, la malla no es mainfold" << endl;
			}
			//Si la arista inversa está en externalEdges
			auto indexInverse = externalAux.find(ed.reversed());
			if (indexInverse != externalAux.end()) {
				//Borrar arista inversa de externalEdges
				externalAux.erase(indexInverse);
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
			k = i;
		// intercambiar pos
		SimpleEdge aux = externalEdges[i];
		externalEdges[i] = externalEdges[k];
		externalEdges[k] = aux;
	}
    //END TODO 2.1

}//void updateEdgeLists()
//////////////////////e
double computeAngle(vec3 &a, vec3 &b, vec3 &c)
{
	vec3 ab = b - a;
	vec3 ac = c - a;

	double m = sqrt((ab*ab)*(ac*ac));
	double cosine = (ab*ac) / m;
	double sine = (ab^ac).module() / m;
	return atan2(sine, cosine);
}
//////////////////////////////7

/// Write an Eigen Matrix to a matlab file
void exportDenseToMatlab (const Eigen::MatrixXd &m, const std::string &filename, const std::string matrixName ="A")
{
    cout << "Export matrix "<< matrixName << " to file: " << filename << std::endl;

    //Open the file as a stream
    ofstream os(filename.c_str());
    if (!os.is_open())
        throw filename + string(": Error creating the file");


    os << "# name: "<<matrixName << std::endl
       << "# type: matrix" << std::endl
       << "# rows: " << m.rows() << std::endl
       << "# columns: " << m.cols() << std::endl;
    Eigen::IOFormat fmt(Eigen::StreamPrecision, Eigen::DontAlignCols, " ", "\n", "", "", "", "\n");
    os << m.format(fmt);

    os.close();
    std::cout << "To import "<< matrixName <<" into matlab use the command: load(\""<<filename<<"\")"<<std::endl;

}//void exportDenseToMatlab (const &MatrixXd m, const std::string &filename)


int main (int argc, char *argv[])
{
    //Set dumpMatrix to true for debugging. Writting matrix to file can take some time
    bool dumpMatrix = false;

    try
    {
        // Set default input mesh filename
		//std::string filename("mallas/16Triangles.off");  //Minimal case test
        std::string filename("mallas/mask2.ply");      //Easy case test
        //std::string filename("mallas/mannequin2.ply"); //Medium case test
        //std::string filename("mallas/laurana50k.ply"); //Really hard for dense matrix
		//std::string filename("mallas/Hexagon2.off");
		//std::string filename("mallas/Lion.ply");
        if (argc > 1)
            filename = std::string(argv[1]);

        ///////////////////////////////////////////////////////////////////////
        //Step 1.
        //Read an input mesh
        SimpleMesh mesh;
        cout << "Loading file " << filename << endl;
        mesh.readFile(filename, false);

        cout << "Num vertex: " << mesh.numVertex() << " Num triangles: " << mesh.numTriangles()
             << " Unreferenced vertex: " << mesh.checkUnreferencedVertex() << endl;

        //Set dumpMatrix to true for debugging. Writting matrix to file can take some time
        dumpMatrix = dumpMatrix || (mesh.numVertex() < 40);

        //Time measure
        high_resolution_clock::time_point clock0 = high_resolution_clock::now();

        ///////////////////////////////////////////////////////////////////////
        //Step 2.
        //Compute edge list and show it. This step should work if you correctly
        //finished the meshBoundary exercise.
        std::vector<SimpleEdge> externalEdges;
        std::vector<SimpleEdge> internalEdges;
        updateEdgeLists(mesh, externalEdges, internalEdges);

        //Count num of internal and external vertex
        size_t numVertex = mesh.coordinates.size();
        size_t numOfExternalVertex = externalEdges.size();
        //size_t numOfInternalVertex = numVertex - numOfExternalVertex;


        //Dump external edges
        cout << numOfExternalVertex << " vertex in external boundary: " << endl;
        if (numOfExternalVertex < 80)
        {
            for (auto &e : externalEdges)
            {
                cout << "[" << e.a << "->" << e.b << "] ";
            }
            cout << endl;
        }

        //Build the list of external vertex
        std::vector<size_t>externalVertex;
        externalVertex.reserve(externalEdges.size());
        for (auto &e : externalEdges)
        {
            externalVertex.push_back(e.a);
        }

        //Optimization: Keep a lookup table to ask if one vertex is external
        std::vector<bool>isExternal(numVertex, false);
        for(size_t i=0; i< numOfExternalVertex; i++)
        {
            isExternal[externalVertex[i]] = true;
        }

        cout << "Done compute Boundary. " << duration<float>(high_resolution_clock::now() - clock0).count() << " seconds" << endl;

        ///////////////////////////////////////////////////////////////////////
        //Step 3.
        //Build Laplace matrix (step 1)
        cout << "Computing Laplacian matrix (dense):" << endl;

        //We will use a dense matrix. If you want to use sparse matrix
        //you have to check http://eigen.tuxfamily.org/dox/group__TutorialSparse.html
        //and declare a sparse matrix instead.
        Eigen::MatrixXd meshMatrix(numVertex, numVertex);

        //TODO 3.1: Build the Laplacian matrix using the weights
        //
        //      /\         .  Lij = cot(α) + cot(β) , if vertex i neighbour of j
        //     /β \        .
        //    /    \       .  Lii = -Sum Lij , diagonal element is the sum of row
        //   /      \      .
        // vi--------vj    .
        //   \      /      .
        //    \    /       .
        //     \α /        .
        //      \/         .
        //
		std::vector<std::vector<int>> neighbour; 
		mesh.computeNeighbours(neighbour);//calcula los vecinos de cada vértice
	
		float solCot = 0;
		float angle;
		
		//recorre la matriz, si el vértice i es vecino del vértice j, Lij = cot(α) + cot(β).
		//Si no, es 0.

		for (int i = 0; i < numVertex; ++i) {
			for (int j = 0; j < numVertex; ++j) {
				//Si el vertice i es vecino del vertice j
				if (std::find(neighbour[i].begin(), neighbour[i].end(), j) != neighbour[i].end())
				{
					//Recorriendo los vecinos del vertice i, miramos si i[k] esta en los vecions del vertice j
					for (int k = 0; k < neighbour[i].size(); ++k) {
						if ((std::find(neighbour[j].begin(), neighbour[j].end(), neighbour[i][k]) != neighbour[j].end()))
						{
							//Calcular ángulo con los 3 puntos del triángulo
							angle = computeAngle(mesh.coordinates[neighbour[i][k]], mesh.coordinates[i], mesh.coordinates[j]);
							solCot += 1 / tan(angle); //guarda la sol, sumaría 2 veces max
						}
					}
					meshMatrix(i, j) = solCot;
					solCot = 0; //reiniciar a 0 para el siguiente vértice/posición de la matriz
				}
				else
					meshMatrix(i, j) = 0;
			}
		}
		
	
		//Calcular los elementos de la diagonal de forma que la suma de los elementos de cada fila
		//valga 0. para ello uso meshMatrix.row(i).sum(), que accede a la fila i y suma sus elementos
		for (int i = 0; i < numVertex; ++i) 
			meshMatrix(i, i) = - meshMatrix.row(i).sum();
			
		
		
	
        //END TODO 3.1

        cout << "Done Laplacian matrix (dense). " << duration<float>(high_resolution_clock::now() - clock0).count() << " seconds" << endl;

        if (dumpMatrix)
            exportDenseToMatlab(meshMatrix, "laplacian.mat", "L");

        ///////////////////////////////////////////////////////////////////////
        //Step 5.1
        //Build system matrix
        cout << "Computing System matrix:" << endl;

        //TODO 3.2: Patch Laplace matrix to generate a valid system of equations
		
		//recorre la matriz, si el vértice i no esta en la frontera lo deja igual
		// si i != j el valor es 0, si i==j el valor es 1.
		for (int i = 0; i < numVertex; ++i) {
			for (int j = 0; j < numVertex; ++j) {
				if (std::find(externalVertex.begin(), externalVertex.end(), i) != externalVertex.end())
				{
					if(i!=j)
						meshMatrix(i, j) = 0;
					else
						meshMatrix(i, j) = 1;
				}
			}
		}
		
		
        //END TODO 3.2

        cout << "Done System matrix. " << duration<float>(high_resolution_clock::now() - clock0).count() << " seconds" << endl;

        if (dumpMatrix)
            exportDenseToMatlab(meshMatrix, "systemMatrix.mat", "A");

        ///////////////////////////////////////////////////////////////////////
        //Step 5.2
        //Build UV values for external vertex, using the boundary of a square.
        cout << "Computing contour conditions:" << endl;

        //Note that this is a matrix of numvertex rows and 2 columns 
        Eigen::MatrixX2d UV_0(numVertex,2);
        UV_0.setZero();
        
        //TODO 3.3: Build a set of valid UV values for external vertex, mapping
        //the vertex to the boundary of a square of side unit.

		double p = 4 / double(numOfExternalVertex); //tamaño de los pasos en los que dividir 
		double u = 0; //valor u
		double v = 0; //Valor v
		int pos = 0; 
		int ct = 0; //índice que recorre externalVertex
	

		//para todos los vertices, se miran los que están en el exterior y se mapea
		for (int i = 0; i < numVertex; ++i) {
			auto id= std::find(externalVertex.begin(), externalVertex.end(), i);
			if ( id != externalVertex.end())       //Si está en el exterior
			{	
				pos = externalVertex[ct]; //Vértice que esta en la posicion ct de externalVertex. 
				//Va asignando valores UV. Dándole valores más bajos a los vértices externos siguiendo el
				//orden en el que están en el vector externalVertex. 
				//Primero mapea en el lado inferior del cuadrado
				if ((u < 1) && (v == 0)) { 
					UV_0(pos, 0) = u;
					UV_0(pos, 1) = v;
					u += p;
					
					//ha llegado al extremo, o casi
					if (u > 1) {
						v += u-1;
						u = 1;
					}
				//Mapeo lado derecho
				}else if ((v < 1) && (u >= 1)) {
					UV_0(pos, 1) = v; 
					UV_0(pos, 0) = u;
					v += p;
					
					if (v > 1) {
						u -= (v-1);
						v = 1;
					}
				//mapeo lado superior
				}else if ((u > 0) && (v >= 1)) {
					UV_0(pos, 0) = u;
					UV_0(pos, 1) = v;
					u -= p;
					
					if (u < 0) {
						v += u;
						u = 0;
					}
				//Mapeo lado izquierdo
				}else if ((v > 0) && (u <= 0)){
					UV_0(pos, 1) = v;
					UV_0(pos, 0) = u;
					v -= p;					
				}				
				++ct; //se pasa al siguiente valor de externalVertex
			}
		}
        //END TODO 3.3
		

        cout << "Done contour conditions. " << duration<float>(high_resolution_clock::now() - clock0).count() << " seconds" << endl;

        if (dumpMatrix)
            exportDenseToMatlab(UV_0, "boundaryUVO.mat", "UV0");

        ///////////////////////////////////////////////////////////////////////
        //Step 6

        //We will use Eigen to solve the matrix from here. 
        cout << "Solving the system using Eigen (dense):" << endl;

        //Solve Ax = b; where A = meshMatrix and b = UV_0
        //This code is valid only for solving a dense matrix. If you want to use sparse matrix
        //you have to check http://eigen.tuxfamily.org/dox/group__TopicSparseSystems.html
        Eigen::PartialPivLU<Eigen::MatrixXd>meshMatrixLU = meshMatrix.lu();
        Eigen::MatrixX2d UV = meshMatrixLU.solve(UV_0);

        cout << "Done solving the system. " << duration<float>(high_resolution_clock::now() - clock0).count() << " seconds" << endl;

        //Dump the computed solution to the system
        if (dumpMatrix)
            exportDenseToMatlab(UV, "solutionUV.mat", "UV");

        ///////////////////////////////////////////////////////////////////////
        //Step 6.4 

        //Create a planar mesh (reverse UV mesh)
        SimpleMesh planarMesh;
        //Use UV values as geometrical coordinates and same triangles that input mesh
        planarMesh.coordinates.resize(numVertex);
        for (size_t i=0; i< numVertex; i++)
        {
            planarMesh.coordinates[i].set(UV(i,0), UV(i,1), 0);
        }
        planarMesh.triangles = mesh.triangles;

        string output_UVMesh1="output_UVMesh1.ply";
        cout << "Saving parameterization mesh to " << output_UVMesh1 << endl;
        //planarMesh.writeFileOBJ(output_UVMesh1);
        planarMesh.writeFilePLY(output_UVMesh1);


        //Create a TextureMesh mesh with UV coordinates
        TextureMesh textureMesh;
        textureMesh.coordinates = mesh.coordinates;
        textureMesh.triangles= mesh.triangles;

        //Set image filename to be used as texture
        textureMesh.textureFile= "UVchecker.jpg";

        //Set UV as texture-per-vertex coordinates
        textureMesh.UV.resize(numVertex);
        for (size_t i=0; i< numVertex; i++)
            textureMesh.UV[i].set(float(UV(i,0)), float(UV(i,1)));

        //Dump textureMesh to file (.obj or .ply)
        string output_UVMesh2="output_UVMesh2.ply";
        cout << "Saving texture mesh to " << output_UVMesh2 << endl;
        //textureMesh.writeFileOBJ(output_UVMesh2);
        textureMesh.writeFilePLY(output_UVMesh2);

        //Visualize the file with an external viewer
#ifdef WIN32
        //string viewcmd = "\"C:\\Program Files (x86)\\VCG\\MeshLab\\meshlab.exe\"";
        string viewcmd = "C:/meshlab/meshlab.exe";
#else
        string viewcmd = "meshlab >/dev/null 2>&1 ";
#endif
        string cmd = viewcmd+" "+output_UVMesh2;
        cout << "Executing external command: " << cmd << endl;
        return system(cmd.c_str());

    }//try
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

