/*
 * meshStatistic.cpp
 *
 * Written by Jose Miguel Espadero <josemiguel.espadero@urjc.es>
 *
 * This code is written as material for the FMF class of the
 * Master Universitario en Informatica Grafica, Juegos y Realidad Virtual.
 * Its purpose is to be didactic and easy to understand, not hard optimized.
 *
 * This file compute some statistics about a mesh and dump then to the console.
 * Also produce an output mesh with degenerate triangles remarked.
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
#include <iomanip>
#include <cstdlib>
#include <SimpleMesh.hpp>
#include <ColorMesh.hpp>

//devuelve el ángulo que forman las aristas a->b y a->c
double computeAngle(vec3 &a, vec3 &b, vec3 &c)
{
	vec3 ab = b - a;
	vec3 ac = c - a;

	double m = sqrt((ab*ab)*(ac*ac));
	double cosine = (ab*ac) / m;
	double sine = (ab^ac).module() / m;
	return atan2(sine, cosine);
}
//devuelve el mínimo de a,b y c
double minimo(double a, double b, double c)
{
	double min;
	if (a < b && a < c)
		min = a;
	if (b < a && b < c)
		min = b;
	if (c < a && c < b)
		min = c;
	return min;
}
//devuelve el máximo de a,b y c
double maximo(double a, double b, double c)
{
	double max;
	if (a > b && a > c)
		max = a;
	if (b > a && b > c)
		max = b;
	if (c > a && c > b)
		max = c;
	return max;
}

//Devuelve si el triangulo es obtuso (para parte opcional)
bool triangleIsObtuse(vec3 &a, vec3 &b, vec3 &c)
{
	//Compute vectors for edges of the triangle
	auto e_ab = b - a;
	auto e_bc = c - b;
	auto e_ca = a - c;

	//Compute cosines of each angle in the triangle.
	double cos_A = -(e_ab * e_ca);
	double cos_B = -(e_bc * e_ab);
	double cos_C = -(e_ca * e_bc);

	//Check directly the cosine (negative for obtuse angles)
	bool obtuse = cos_A < 0 || cos_B < 0 || cos_C < 0;
	return obtuse;
}
//Devuelve si el angulo en el vertice a es obtuso (para parte opcional)
bool angleIsObtuse(vec3 &a, vec3 &b, vec3 &c)
{
	//Compute vectors for edges
	auto e_ab = b - a;
	auto e_ca = a - c;

	//Compute cosine
	double cos_A = -(e_ab * e_ca);


	//Check directly the cosine (negative for obtuse angles)
	bool obtuse = cos_A < 0;
	return obtuse;
}
int main (int argc, char *argv[])
{
    try
    {
        //Set default input mesh filename
        std::string filename("mallas/mannequin2.ply");
		//std::string filename("mallas/mask2.off");
		//std::string filename("mallas/16Triangles.off");
		//std::string filename("mallas/bunny.ply");
		//std::string filename("mallas/AthenaBust.ply");
        if (argc >1)
            filename = std::string(argv[1]);

        //Set default degenerate triangle shapeFactor threshold
        double shapeFactorTh = 1/(4*sqrt(3));
        if (argc > 2)
            shapeFactorTh = strtod(argv[2], nullptr);

        ///////////////////////////////////////////////////////////////////////
        //Read a mesh from given filename
        SimpleMesh mesh;
        cout << "Loading file " << filename << endl;
        mesh.readFile(filename);

        cout << "\nNum vertex: " << mesh.numVertex() << " Num triangles: " << mesh.numTriangles()
             << " Unreferenced vertex: " << mesh.checkUnreferencedVertex() << endl;


        //Compute minimum, maximum and average area per triangle, and total area for the mesh
        double minArea, maxArea, averageArea=0, totalArea;

        //Compute minimum and maximum angle for the mesh
        double minAngle, maxAngle;

        //Compute minimum, maximum and average edge length
        double minEdgeLen, maxEdgeLen, averageEdgeLen;

        //Compute minimum, maximum and average shape factor
        double minShapeFactor, maxShapeFactor, averageShapeFactor;

        //////////////////////////////////////////////////////////////////////////////////
        //TODO 1.1:
        // Compute the values for minArea, maxArea, averageArea, totalArea
		double areaActual;
		minArea = INFINITY;
		maxArea = 0;
		averageArea = 0;
		totalArea = 0;
		for (int t = 0; t < mesh.numTriangles(); ++t) {
			areaActual = mesh.triangleArea(mesh.triangles[t]);
			totalArea += areaActual;
			if (areaActual < minArea)
				minArea = areaActual;
			if (areaActual > maxArea)
				maxArea = areaActual;
		}
		averageArea = totalArea / mesh.numTriangles(); 


        // Compute the values for minAngle, maxAngle;
		minAngle = INFINITY;
		maxAngle = 0;
		double actualAngleA;
		double actualAngleB;
		double actualAngleC;
		for (int t = 0; t < mesh.numTriangles(); ++t) {
			actualAngleA = computeAngle(mesh.coordinates[mesh.triangles[t].a], mesh.coordinates[mesh.triangles[t].b], mesh.coordinates[mesh.triangles[t].c]);
			actualAngleB = computeAngle(mesh.coordinates[mesh.triangles[t].b], mesh.coordinates[mesh.triangles[t].c], mesh.coordinates[mesh.triangles[t].a]);
			actualAngleC = computeAngle(mesh.coordinates[mesh.triangles[t].c], mesh.coordinates[mesh.triangles[t].a], mesh.coordinates[mesh.triangles[t].b]);
			
			if (minimo(actualAngleA, actualAngleB, actualAngleC) < minAngle)
				minAngle = minimo(actualAngleA, actualAngleB, actualAngleC);
			if (maximo(actualAngleA, actualAngleB, actualAngleC) > maxAngle)
				maxAngle = maximo(actualAngleA, actualAngleB, actualAngleC);
			//cout << "Minimo " << minimo(actualAngleA, actualAngleB, actualAngleC) << ", min: " << min(actualAngleA, actualAngleB, actualAngleC) << endl;

		}


        // Compute the values for minEdgeLen, maxEdgeLen, averageEdgeLen;
		double edgeTotal=0;
		minEdgeLen = INFINITY;
		maxEdgeLen = 0;
		averageEdgeLen = 0;
		double minEdgeActual;
		double maxEdgeActual;
		int contEdges = 0;
		for (int t = 0; t < mesh.numTriangles(); ++t) {
			for (int e = 0; e < 3; ++e) {
				++contEdges;
				edgeTotal+= mesh.edgeLength(mesh.triangles[t].a,mesh.triangles[t].b)+ mesh.edgeLength(mesh.triangles[t].b, mesh.triangles[t].c)+ mesh.edgeLength(mesh.triangles[t].c, mesh.triangles[t].a);
				minEdgeActual = mesh.triangleMinEdge(mesh.triangles[t]);
				maxEdgeActual = mesh.triangleMaxEdge(mesh.triangles[t]);
				if (minEdgeActual < minEdgeLen)
					minEdgeLen = minEdgeActual;
				if (maxEdgeActual > maxEdgeLen)
					maxEdgeLen = maxEdgeActual;
			}
		}
		averageEdgeLen = edgeTotal / (3*contEdges);

        // Compute the values for minShapeFactor, maxShapeFactor, averageShapeFactor;
		minShapeFactor = 1;
		maxShapeFactor = 0;
		averageShapeFactor = 0;
		double totalSF = 0;
		double actualShapeFactor;
		for (int t = 0; t < mesh.numTriangles(); ++t) {
			actualShapeFactor = mesh.triangleShapeFactor(mesh.triangles[t]);
			totalSF += actualShapeFactor;
			if (actualShapeFactor < minShapeFactor)
				minShapeFactor = actualShapeFactor;
			if (actualShapeFactor > maxShapeFactor)
				maxShapeFactor = actualShapeFactor;
		}
		averageShapeFactor= totalSF / mesh.numTriangles();
        //END TODO 1.1
		
        // Dump statistics to console output
        cout << std::fixed << std::setprecision(4) <<
                "\nArea    min: " << minArea <<
                " max: " << maxArea  <<
                " average: " << averageArea  <<
                " total: " << totalArea <<
                "\nAngle   min: " << minAngle <<
                " max: " << maxAngle  <<
                "\nEdgeLen min: " << minEdgeLen <<
                " max: " << maxEdgeLen  <<
                " average: " << averageEdgeLen  <<
                "\nShapeF  min: " << minShapeFactor <<
                " max: " << maxShapeFactor  <<
                " average: " << averageShapeFactor  << endl;

        //////////////////////////////////////////////////////////////////////////////////
        //TODO OPTATIVE 1:
        //Compute Vertex area for each vertex and store in vertexAreas vector
        //Compute minimun and maximun vertex area in minVertexArea, maxVertexArea
        std::vector<double>vertexAreas;
        double sumVertexAreas = 0;
        double minVertexArea = INFINITY;
        double maxVertexArea = 0;
		double actVertexArea = 0;
		
		std::vector<std::vector<int>> neighbour;
		mesh.computeNeighbours(neighbour);//calcula los vecinos de cada vértice

		for (int i = 0; i < mesh.numVertex(); ++i) {
			for (int j = 0; j < mesh.numTriangles(); ++j) {
				//Si T es obtuso:
				if (mesh.triangles[j].a == i|| mesh.triangles[j].b == i || mesh.triangles[j].c == i) {				
						auto vertA =( mesh.coordinates[mesh.triangles[j].a]);
						auto vertB = mesh.coordinates[mesh.triangles[j].b];
						auto vertC = mesh.coordinates[mesh.triangles[j].c];
					 if (mesh.triangles[j].b == i) {
						auto vertA = mesh.coordinates[mesh.triangles[j].b];
						auto vertB = mesh.coordinates[mesh.triangles[j].c];
						auto vertC = mesh.coordinates[mesh.triangles[j].a];
					}else if (mesh.triangles[j].c == i) {
						auto vertA = mesh.coordinates[mesh.triangles[j].c];
						auto vertB = mesh.coordinates[mesh.triangles[j].a];
						auto vertC = mesh.coordinates[mesh.triangles[j].b];
					}
					if (triangleIsObtuse(vertA, vertB, vertC)) {
						if (angleIsObtuse(vertA, vertB, vertC)) {
							actVertexArea += (mesh.triangleArea(mesh.triangles[j]) / 2.0);
						}
						else {
							actVertexArea += (mesh.triangleArea(mesh.triangles[j]) / 4.0);
						}
					}
					else {
						//veronoi es seguro
						actVertexArea += (1.0 / 8.0)* (pow((vertB - vertA).module(), 2.0)* (1.0 / tan(computeAngle(vertC, vertA, vertB)))
							+ pow((vertC - vertA).module(), 2.0)* (1.0 / tan(computeAngle(vertB, vertC, vertA))));
					}
				}
			}
			
			vertexAreas.push_back(actVertexArea);
			sumVertexAreas += actVertexArea;
			if (actVertexArea < minVertexArea)
				minVertexArea = actVertexArea;
			if (actVertexArea > maxVertexArea)
				maxVertexArea = actVertexArea;
			actVertexArea = 0;
		}


        //END TODO OPTATIVE 1

        //Check values for vertex areas and compute sum of areas:
        if (vertexAreas.size() != mesh.numVertex() )
        {
          cout << "VertexAreas OPTATIVE PART NOT DONE" << endl;
        }
        else
        {
          cout << "VertexA min: " << minVertexArea <<
                  " max: " << maxVertexArea  <<
                  " average: " << sumVertexAreas / mesh.numVertex()  <<
                  " total: " << sumVertexAreas << endl;
        }


        //////////////////////////////////////////////////////////////////////////////////
        //TODO 1.2:
        //Write triangles with ShapeFactor (radius / minEdge) greater than shapeFactorTh
        //Use messages formated as: "Triangle nnn has ShapeFactor xxx"
		cout << endl;
		for (int t = 0; t < mesh.numTriangles(); ++t) {
			actualShapeFactor = mesh.triangleShapeFactor(mesh.triangles[t]);
			if (actualShapeFactor < shapeFactorTh)
				cout << "Triangle " << t << " has ShapeFactor " << actualShapeFactor << endl;
		}


        //END TODO 1.2

        //////////////////////////////////////////////////////////////////////////////////
        //TODO 1.3:
        //Create a colorMesh where faces with ShapeFactor greater than shapeFactorTh
        //have their vertex colored in red. Save it to file named output_statistic.ply
        //and visualize it with meshlab or another external viewer.
		//FALTA
		
		ColorMesh colorMesh;
		//Prepare a copy of the input mesh on a ColorMesh
		colorMesh.coordinates = mesh.coordinates;
		colorMesh.triangles = mesh.triangles;
		colorMesh.colors.resize(mesh.numVertex());
		//recorre los triangulos, y pinta de rojo los que tengan factor de forma inferior a shapeFactorTh
		for (int t = 0; t < colorMesh.numTriangles(); ++t) {
			actualShapeFactor = colorMesh.triangleShapeFactor(colorMesh.triangles[t]);
			if (actualShapeFactor < shapeFactorTh) {
				//vec3 re = colorMesh.coordinates[colorMesh.triangles[t].a];
				colorMesh.colors[mesh.triangles[t].a].set(1.0f, 0.0f, 0.0f);
				colorMesh.colors[mesh.triangles[t].b].set(1.0f, 0.0f, 0.0f);
				colorMesh.colors[mesh.triangles[t].c].set(1.0f, 0.0f, 0.0f);
			}
				
		}
	

		//Save result to a file in .ply format
		string outputFilename = "output_meshStatistic.ply";
		cout << "Saving output to " << outputFilename << endl;
		colorMesh.writeFilePLY(outputFilename);

		//Visualize the .ply file with an external viewer
		string viewcmd = "C:/meshlab/meshlab.exe";
		string cmd = viewcmd + " " + outputFilename;
		cout << "Executing external command: " << cmd << endl;
		return system(cmd.c_str());
        //END TODO 1.3
		
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

