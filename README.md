# A Sagemath tool for abelian surface simplicity
Source code for the Sagemath based software to find abelian surface simplicity
Program for searching for abelian submanifolds of order 2 associated with a product of elliptic curves.
To begin simply download 'surface2.sage' and use load('surface2.sage') to load the program to your notebook. (If surface2.safe is not on the same folder as your notebook you must specify the path to it.)
The code is built around the Surface2 class.
To begin, a Surface2(d2,t1,t2,t3) object is created using the matrix parameters that define the system of 3 equations.

S = Surface2(d2,t1,t2,t3)

Using this object you can use all the functions provided here, from finding solutions to the sistem to check if your input is valid.
The principal method of this class is 'resolution', used to search for solutions to the system of equations. The search algorythms turn the equations into an ellipsoid surface and evaluat it in a discrete grid to take and compare samples of the it.

S.resolution(order, method = "ellipsoid_looping", point_limit = -1)

The input parameters are as follows:

order: determines how thin is the grid from which samples are taken in the space of the surface. A sample will be taken each 10^-order. A greater value will evaluate more samples ang give more results or better chances to find a result at the cost of increased execution time.

method: Choses one of the two methods avilable for the seach. The default methos is 'ellipsoid_looping', it takes samples around the areas in which solutions can be found and only searches around them. It's the fastest method, but there is a chance it may become unstable or not find a certain solution or that it finds the same solution twice in certain edge cases.
The second method is called 'full_scan'. This method makes a full sweep in the area that bounds the ellipsoid with the solutions. If there are solutions with a certain given order this method will find them, with the disadvantage that it is slower and scales with notation O(n^3) with n being 'order'.

point_limit: This is a value to limit the ammount of samples taken by the program so that the user may controll the execution time. If no value is given, the value 10000^order is used.

The expected output is a queue from the 'collections' library with all the 6 value tuples that solve the system. There will also be printed a 3D graph of the ellipsoid surface with the solutions, the 6 parameter tuples that solve the system and the ammount of solutions found.

Another usefull method is 'ascertainment'

S.ascertainment(a12,a13,a14,a23,a24,a34)

This method takes the 6 values of a tuple that is candidate for solution to the system, evaluates them in its 3 equations and returns True if the tuple is solution and false otherwise.
