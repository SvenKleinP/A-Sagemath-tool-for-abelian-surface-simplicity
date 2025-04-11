from operator import itemgetter
from collections import deque

class Surface2:
    def __init__(self,d2,t1,t2,t3):
        self.d2 = d2
        self.t1r = t1.real()
        self.t1i = t1.imag()
        self.t2r = t2.real()
        self.t2i = t2.imag()
        self.t3r = t3.real()
        self.t3i = t3.imag()
        self.aa = (-self.t2i*self.t2i*self.t3i + self.t2r*self.t2r*self.t3i + self.t1i*self.t3i*self.t3i 
                   - 2*self.t2i*self.t2r*self.t3r + self.t1i*self.t3r*self.t3r)/(self.t3i*d2)
        self.bb = d2
        self.cc = d2*self.t1i/self.t3i
        self.dd = -2*d2*self.t2i/self.t3i
        self.ee = 2*(self.t2i*self.t2r - self.t1i*self.t3r)/self.t3i
        self.ff = 2*(self.t2i*self.t3r - self.t2r*self.t3i)/self.t3i
        self.gg = (self.t2i*self.t3r - self.t2r*self.t3i)/self.t3i
        self.hh = d2
        self.ii = -d2*self.t2i/self.t3i
        
#@brief checks if the surface matrix is positive
#@return True if the surface matrix is positive. False it is not
    def selfCheck(self):
        if(self.aa<=0):
            return False
        if(self.aa*self.bb-pow(self.ff/2,2)<=0):
            return False
        if(self.aa*self.bb*self.cc + self.dd*self.ee*self.ff/4 - self.bb*pow(self.ee/2,2) - self.aa*pow(self.dd/2,2) - self.cc*pow(self.ff/2,2) <= 0):
            return False
        if(self.t1i<=0):
            return False
        if(self.t1i*self.t3i-pow(self.t2i,2)<=0):
            return False
        return True
    
#@brief method that calculates the value of the surface function for a given 3 variable tuple
#@param a12, a13 and a14 are the values that are being evaluated in the equation
#@return the value of the surface function for a given a12, a13 and a14 input
    def calculate(self,a12,a13,a14):
        return (a12*a12*self.aa + a13*a13*self.bb + a14*a14*self.cc + a13*a14*self.dd + a12*a14*self.ee + a12*a13*self.ff + a12*self.gg + a13*self.hh + a14*self.ii)
    
#@brief: method that returns the surface equations with variables 'var' a12, a13 y a14
#@return: the surface equationde eq(a12,a13,a14)
    def equation(self):
        a12,a13,a14 = var('a12,a13,a14')
        return a12*a12*self.aa + a13*a13*self.bb + a14*a14*self.cc + a13*a14*self.dd + a12*a14*self.ee + a12*a13*self.ff + a12*self.gg + a13*self.hh + a14*self.ii
    
#@brief method to find the maximum and minimum in each axle of the figure
#functions finding the gradient of the figure that is in the direction of each axle.
#@return 3*2*3array with the coordinates of the maximums and minimums of a12, a13 and a14 [[a12min,a12max],[a13min,a13max],[a14min,a14max]]
    def bounds(self):
        var('x y z')
        minmaxx = solve((x^2*self.aa + y^2*self.bb + z^2*self.cc + y*z*self.dd + x*z*self.ee + x*y*self.ff + x*self.gg + y*self.hh + z*self.ii==0,
                        2*y*self.bb + z*self.dd + x*self.ff + self.hh==0,
                        2*z*self.cc + y*self.dd + x*self.ee + self.ii==0),x,y,z)
        minmaxy = solve((2*x*self.aa + z*self.ee + y*self.ff + self.gg==0,
                         x^2*self.aa + y^2*self.bb + z^2*self.cc + y*z*self.dd + x*z*self.ee + x*y*self.ff + x*self.gg + y*self.hh + z*self.ii==0,
                        2*z*self.cc + y*self.dd + x*self.ee + self.ii==0),x,y,z)
        minmaxz = solve((2*x*self.aa + z*self.ee + y*self.ff + self.gg==0,
                        2*y*self.bb + z*self.dd + x*self.ff + self.hh==0,
                        x^2*self.aa + y^2*self.bb + z^2*self.cc + y*z*self.dd + x*z*self.ee + x*y*self.ff + x*self.gg + y*self.hh + z*self.ii==0),x,y,z)
#this sorts the bound of each axle to be the minumum and maximum, in that order
        limitesx = sorted([[minmaxx[0][0].rhs(),minmaxx[0][1].rhs(),minmaxx[0][2].rhs()],
                [minmaxx[1][0].rhs(),minmaxx[1][1].rhs(),minmaxx[1][2].rhs()]],key = itemgetter(0))
        limitesy = sorted([[minmaxy[0][0].rhs(),minmaxy[0][1].rhs(),minmaxy[0][2].rhs()],
                [minmaxy[1][0].rhs(),minmaxy[1][1].rhs(),minmaxy[1][2].rhs()]],key = itemgetter(1))
        limitesz = sorted([[minmaxz[0][0].rhs(),minmaxz[0][1].rhs(),minmaxz[0][2].rhs()],
                [minmaxz[1][0].rhs(),minmaxz[1][1].rhs(),minmaxz[1][2].rhs()]],key = itemgetter(2))
        return[limitesx,limitesy,limitesz]
    
#@brief method to undo the simplification process. Returns the 6 variable array
#@param a12, a13 and a14 are variables for a given sub manyfold
#@return an array with the variables [a12,a13,a14,a23,a24,a34]
    def reverse(self,a12,a13,a14):
        a24 = -self.d2 -a13*self.d2
        a23 = (a14*self.d2*self.t1i - self.d2*self.t2i - 2*a13*self.d2*self.t2i + 2*a12*self.t2i*self.t2r - a12*self.t1r*self.t3i - 
 a12*self.t1i*self.t3r)/self.t3i
        a34 = (a14*self.d2*self.t1r*self.t3i - a12*pow(self.t2i,2)*self.t3i - self.d2*self.t2r*self.t3i - 2*a13*self.d2*self.t2r*self.t3i + 
 a12*pow(self.t2r,2)*self.t3i + a12*self.t1i*pow(self.t3i,2) - a14*self.d2*self.t1i*self.t3r + self.d2*self.t2i*self.t3r + 
 2*a13*self.d2*self.t2i*self.t3r - 2*a12*self.t2i*self.t2r*self.t3r + a12*self.t1i*pow(self.t3r,2))/(self.d2*self.t3i)
        return [a12,a13,a14,a23,a24,a34]
    
#@brief Plots in 3D the sub manyfold in the coordinates a12, a13 y a14
#@param 3*2 array of the superior and inferiorlimits of a12,a13 y a14
#[[a12min,a12max],[a13min,a13max],[a14min,a14max]]
#@return a 3D plot of the ellipsoid.
    def plt(self,limites):
        var('a12 a13 a14')
        return implicit_plot3d(self.aa*a12^2 + self.bb*a13^2 + self.cc*a14^2 + a13*a14*self.dd + a12*a14*self.ee + a12*a13*self.ff + a12*self.gg + a13*self.hh + a14*self.ii==0,
                               (a12,limites[0][0],limites[0][1]),(a13,limites[1][0],limites[1][1]),(a14,limites[2][0],limites[2][1]))

#@brief Evaluates the 3 original equations with the input variable values for the given sub manyfold
#@param a12, a13, a14, a23, a24 and a34 variables for the tuple to evaluate. functions better inputing them in fraction form
#@return returns True if the 3 equations are exactly equal to 0. False in any other case.
    def ascertainment(self,a12,a13,a14,a23,a24,a34):
        a = SR(a23)
        b = SR(a24)
        c = SR(a34)
        eq1 = self.d2*a13+a24+self.d2
        eq2 = (((self.t1r+I*self.t1i)*(self.t3r+I*self.t3i)-(self.t2r+I*self.t2i)^2)*a12
           -self.d2*(self.t1r+I*self.t1i)*a14 + self.d2*(self.t2r+I*self.t2i)*a13
          -(self.t2r+I*self.t2i)*a24 + (self.t3r+I*self.t3i)*a23 + self.d2*a34)
        eq3 = a14*a23-a13*a24+a12*a34
        #the 3 last variables have to be tested to check that they do not contain irrational implicit expressions
        if(not a.is_numeric()or not b.is_numeric() or not c.is_numeric()):
            return False
        if(not eq1.is_zero()):
            return False
        if(not eq2.is_zero()):
            return False
        if(not eq3.is_zero()):
            return False
        return True

#@brief Generates the surface equation, bounds of the solution area and evaluates the surface to find solutions.
#@param order the size of the discretization that will subdivide the space of the system. it is subdivided in spaces of 10^-n
#@param method the method used for the seatch. "ellipsoid_looping" is the default method.
#The "full_scan" method can be selected to test all pints in the bounds of the system.
#@param point_limit sets a maximum number of tests before the search is concluded. If no limit is setm 10000^order is used.
#@return 0 upon finishing. Shows the sub manyfold plot and prints all the obtained results. Returns -1 in case of an error.
    def resolution(self, order, method = "ellipsoid_looping", point_limit = -1):
        if(not self.selfCheck()):
            print("Error in the curve parameters.")
            return -1
        if(order<1| (not isinstance(order,sage.rings.integer.Integer))):
            print("Error in the input parameters.")
            return -1
        spacing = 1/pow(10,order)
        nro_soluciones = 0
        soluciones = deque()
        lim = self.bounds()#the 3*2*3 limt array is remade into a 3*2 array with the maximums and minimums of each variable
        #[[min a12, max a12],[min a13, max a13],[min a14, max a14]]
        limites = [[lim[0][0][0],lim[0][1][0]],[lim[1][0][1],lim[1][1][1]],[lim[2][0][2],lim[2][1][2]]]
        show(self.plt(limites))
        limites = [[QQ(round(limites[0][0],order)),QQ(round(limites[0][1],order))],
                   [QQ(round(limites[1][0],order)),QQ(round(limites[1][1],order))],
                   [QQ(round(limites[2][0],order)),QQ(round(limites[2][1],order))]]
        if(method=="full_scan"):
            x = limites[0][0]
            print('|{:^7}|{:^7}|{:^7}|{:^7}|{:^7}|{:^7}|'.format("a12","a13","a14","a23","a24","a34"))
            while x <= limites[0][1]:
                y = limites[1][0]
                while y <= limites[1][1]:
                    z = limites[2][0]
                    while z <= limites[2][1]:
                        if(self.calculate(x,y,z).is_zero()):
                            answer = self.reverse(x,y,z)
                            if(self.ascertainment(x,y,z,answer[3],answer[4],answer[5])):
                                nro_soluciones+=1
                                soluciones.append([answer[0],answer[1],answer[2],answer[3],answer[4],answer[5]])
                                print('|{:^3}/{:^3}|{:^3}/{:^3}|{:^3}/{:^3}|{:^3}/{:^3}|{:^3}/{:^3}|{:^3}/{:^3}|'.format(x.numerator(),x.denominator(),y.numerator(),y.denominator(),z.numerator(),z.denominator(),QQ(answer[3]).numerator(),QQ(answer[3]).denominator(),QQ(answer[4]).numerator(),QQ(answer[4]).denominator(),QQ(answer[5]).numerator(),QQ(answer[5]).denominator()))
                        z=z+spacing
                    y=y+spacing
                x=x+spacing
            print(nro_soluciones)
            return soluciones
        if(method == "ellipsoid_looping"):
            if(point_limit < 0):
                point_limit = pow(10000,order)
            x = QQ(round(lim[0][0][0],order))
            if(abs(x)>abs(lim[0][0][0])):
                x = x + spacing
            print('|{:^7}|{:^7}|{:^7}|{:^7}|{:^7}|{:^7}|'.format("a12","a13","a14","a23","a24","a34"))
            p1=((QQ(round(lim[0][0][1],order)),QQ(round(lim[0][0][2],order))))#samples are taken in p1 and p2 comparing their values
            p2=((QQ(round(lim[0][0][1],order)),QQ(round(lim[0][0][2],order))))#p1 goes in front of
            n=0
            #this loop moves the point 1 to the border of the circle. after crossingit it goes back to the interior point closer to the border
            while(self.calculate(x,p1[0],p1[1]).abs()<=self.calculate(x,p2[0],p2[1]).abs() and n<point_limit):
                p2=p1
                temp = list(p1)
                temp[0] = temp[0]+spacing
                temp[1] = temp[1]+spacing
                p1=tuple(temp)
                n=n+1
            p1 = p2
            #once that the border is found in the a13+ way, the point 1 is displaced in the a14- way untill the ellipsoid is crossed
            #this is to make sure that the starting point is in the quadrant 1 a13 positive and a14 positive just outside the ellipsoid
            while(self.calculate(x,p1[0],p1[1]).abs()<=self.calculate(x,p2[0],p2[1]).abs() and n<point_limit):
                p2=p1
                temp = list(p1)
                temp[1] = temp[1]-spacing
                p1=tuple(temp)
                n=n+1
            p1=p2
        
            while(x<limites[0][1][0] and n<point_limit):
                k = 0#k marks which quadrant is the current search in.
                direccion = True#determines the variable that is being used to advance, either a13 o a14. when it changes there is a change in direction.
                temp = list(p1)
                temp[0] = temp[0]+spacing
                p1=tuple(temp)
                while(self.calculate(x,p1[0],p1[1]).abs()<=self.calculate(x,p2[0],p2[1]).abs() and self.calculate(x,p1[0],p1[1])*self.calculate(x,p2[0],p2[1])>0):
                    p2=p1
                    temp = list(p1)
                    temp[0] = temp[0]+spacing
                    p1=tuple(temp)
                    n=n+1
                while(k<4 and n<point_limit):#The start of the search is in the border of the first and fourth quadrant and runs counter clockwise
                    p2=p1
                    temp = list(p1)
                    if(direccion):
                        if(k==0):
                            temp[0] = temp[0]-spacing
                        elif(k==1):
                            temp[1] = temp[1]-spacing
                        elif(k==2):
                            temp[0] = temp[0]+spacing
                        elif(k==3):
                            temp[1] = temp[1]+spacing
                    else:
                        if(k==0):
                            temp[1] = temp[1]+spacing
                        elif(k==1):
                            temp[0] = temp[0]-spacing
                        elif(k==2):
                            temp[1] = temp[1]-spacing
                        elif(k==3):
                            temp[0] = temp[0]+spacing
                    p1=tuple(temp)
                    n=n+1
                    if(self.calculate(x,p1[0],p1[1])*self.calculate(x,p2[0],p2[1])<=0):#the frontier of the ellipsoid was crossed
                        direccion = not direccion
                        if(self.calculate(x,p1[0],p1[1]).is_zero()):
                            answer = self.reverse(x,p1[0],p1[1])
                            if(self.ascertainment(x,p1[0],p1[1],answer[3],answer[4],answer[5])):
                                nro_soluciones +=1
                                soluciones.append([answer[0],answer[1],answer[2],answer[3],answer[4],answer[5]])
                                print('|{:^3}/{:^3}|{:^3}/{:^3}|{:^3}/{:^3}|{:^3}/{:^3}|{:^3}/{:^3}|{:^3}/{:^3}|'.format(x.numerator(),x.denominator(),p1[0].numerator(),p1[0].denominator(),p1[1].numerator(),p1[1].denominator(),QQ(answer[3]).numerator(),QQ(answer[3]).denominator(),QQ(answer[4]).numerator(),QQ(answer[4]).denominator(),QQ(answer[5]).numerator(),QQ(answer[5]).denominator()))
                    elif(self.calculate(x,p1[0],p1[1]).abs()>(self.calculate(x,p2[0],p2[1]).abs())):#if the point 1 is bigger than 2, then the systen is going away from the ellipse and a change of quadrant is needed
                        #the following checks are used to evaluate the points in the diagonals of the point when changing quadrants
                        #sometimes this methods skips them as solutions otherwise.
                        if(k==0 and self.calculate(x,p1[0]+spacing,p1[1]+spacing).is_zero()):
                            answer = self.reverse(x,p1[0]+spacing,p1[1]+spacing)
                            if(self.ascertainment(answer[0],answer[1],answer[2],answer[3],answer[4],answer[5])):
                                nro_soluciones +=1
                                soluciones.append([answer[0],answer[1],answer[2],answer[3],answer[4],answer[5]])
                                print('|{:^3}/{:^3}|{:^3}/{:^3}|{:^3}/{:^3}|{:^3}/{:^3}|{:^3}/{:^3}|{:^3}/{:^3}|'.format(answer[0].numerator(),answer[0].denominator(),answer[1].numerator(),answer[1].denominator(),answer[2].numerator(),answer[2].denominator(),QQ(answer[3]).numerator(),QQ(answer[3]).denominator(),QQ(answer[4]).numerator(),QQ(answer[4]).denominator(),QQ(answer[5]).numerator(),QQ(answer[5]).denominator()))
                        elif(k==1 and self.calculate(x,p1[0]+spacing,p1[1]-spacing).is_zero()):
                            answer = self.reverse(x,p1[0]+spacing,p1[1]-spacing)
                            if(self.ascertainment(answer[0],answer[1],answer[2],answer[3],answer[4],answer[5])):
                                nro_soluciones +=1
                                soluciones.append([answer[0],answer[1],answer[2],answer[3],answer[4],answer[5]])
                                print('|{:^3}/{:^3}|{:^3}/{:^3}|{:^3}/{:^3}|{:^3}/{:^3}|{:^3}/{:^3}|{:^3}/{:^3}|'.format(answer[0].numerator(),answer[0].denominator(),answer[1].numerator(),answer[1].denominator(),answer[2].numerator(),answer[2].denominator(),QQ(answer[3]).numerator(),QQ(answer[3]).denominator(),QQ(answer[4]).numerator(),QQ(answer[4]).denominator(),QQ(answer[5]).numerator(),QQ(answer[5]).denominator()))
                        elif(k==2 and self.calculate(x,p1[0]-spacing,p1[1]-spacing).is_zero()):
                            answer = self.reverse(x,p1[0]-spacing,p1[1]-spacing)
                            if(self.ascertainment(answer[0],answer[1],answer[2],answer[3],answer[4],answer[5])):
                                nro_soluciones +=1
                                soluciones.append([answer[0],answer[1],answer[2],answer[3],answer[4],answer[5]])
                                print('|{:^3}/{:^3}|{:^3}/{:^3}|{:^3}/{:^3}|{:^3}/{:^3}|{:^3}/{:^3}|{:^3}/{:^3}|'.format(answer[0].numerator(),answer[0].denominator(),answer[1].numerator(),answer[1].denominator(),answer[2].numerator(),answer[2].denominator(),QQ(answer[3]).numerator(),QQ(answer[3]).denominator(),QQ(answer[4]).numerator(),QQ(answer[4]).denominator(),QQ(answer[5]).numerator(),QQ(answer[5]).denominator()))
                        elif(k==3 and self.calculate(x,p1[0]-spacing,p1[1]+spacing).is_zero()):
                            answer = self.reverse(x,p1[0]-spacing,p1[1]+spacing)
                            if(self.ascertainment(answer[0],answer[1],answer[2],answer[3],answer[4],answer[5])):
                                nro_soluciones +=1
                                soluciones.append([answer[0],answer[1],answer[2],answer[3],answer[4],answer[5]])
                                print('|{:^3}/{:^3}|{:^3}/{:^3}|{:^3}/{:^3}|{:^3}/{:^3}|{:^3}/{:^3}|{:^3}/{:^3}|'.format(answer[0].numerator(),answer[0].denominator(),answer[1].numerator(),answer[1].denominator(),answer[2].numerator(),answer[2].denominator(),QQ(answer[3]).numerator(),QQ(answer[3]).denominator(),QQ(answer[4]).numerator(),QQ(answer[4]).denominator(),QQ(answer[5]).numerator(),QQ(answer[5]).denominator()))
                        k+=1
                x = x+spacing
            print(nro_soluciones, ' solutions.')
            return soluciones
        else:
            print("Error selecting method.")
            return -1
