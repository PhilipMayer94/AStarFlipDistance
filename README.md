# FlipDistAstar
The Code for the paper "Engineering an Astar Algorithm for the Flip Distance of Plane Triangulations".

# Dependencies
C++:

	- CGAL
	- Gurobi 9.5 (If you want to use gurobi 10.0 you may have to change FindGurobi.cmake)
 	- Cgal 5.6

Python:

 	- pandas
	- seaborn
 	- matplotlib



# Our result files
The data is given with our complete result files. If you want to repeat our experiments you first need to execute purge.py.



# Building our code 
```bash
git clone {link here}
```
Then move into the project folder and:
```bash
mkdir build
cd build
cmake ..
make 
```

Then (if you have purged the result files) you can either execute the given shell scripts, or call the binary with four or five arguments:
```bash
./A_star_for_flipdistance_experiments algorithm instances size order 
./A_star_for_flipdistance_experiments algorithm instances size order path
```
For the decomposition experiments use the second binary:
```bash
./A_star_for_flipdistance_decomposition decomposition_runtime_experiments
./A_star_for_flipdistance_decomposition decomposition_runtime_experiments path
```
with:

	1. algorithm = eppstein heuristic combined decomposition ilp bfs heuristic
	2. instances = sealevel random
	3. size = 10 15 20 25 30 for random and 25 30 35 41 for sealevel
	4. order = 2,5, 100 (for unconstrained) for sealevel and always 1 for random
	5. path = path to the folder that contains "random_experiments_paper" and "sealevel_experiments_paper" and "random_hod_experiments_paper"

If you do not specify the path it's assumed that the folder "data" is on the same level as build.





If you want to run the code on your own instance use the third binary:
```bash
./A_star_for_flipdistance vertex_path triangulation_1_path triangulation_2_path algorithm output_path output_file_name
```
with:

	1. vertex_path = path to a list of vertices (every row is corresponds to a vertex) 
	2. triangulation_1_path = path to a file of a triangulation (needs to be compatible  with the vertexlist every row is given by three integers corresponding to a triangle on the vertices)
	3. triangulation_2_path = path to a file of a triangulation (needs to be compatible  with the vertexlist every row is given by three integers corresponding to a triangle on the vertices)
	4. algorithm = eppstein heuristic combined decomposition ilp bfs heuristic
	5. output_path = path to the folder that contains should contain the result json
 	6. output_file_name = name of the resultfile


