package primary.methods;
/**
 * @author heynderickxd
 * @author vannc
 */
import java.util.Scanner;
public class PrimaryMethods {
   /**set of methods for dealing with kinematic equations
    * Note that the methods can be used for both linear and angular kinematics
    * problems. 
    * Simply substitute like this:
    * displacement = theta
    * initial velocity = sigma_o (w thing?)
    * final velocity = sigma_f 
    * acceleration = alfa
    */ 
   //methods for first kinematic equation (x=1/2*a*t^2 + v_ot)
   // solving for acceleration 
    public double kinematics1 (double displacement, double time, 
                            double initialVelocity) {
    double acceleration;
    acceleration = (displacement - (initialVelocity * time)) / Math.pow(time,2.0);
    return acceleration;
    }
   // solving for displacement
    public double kinematics2 (double acceleration, double time, 
                                double initialVelocity) {
    double displacement;
    displacement = (0.5 * acceleration * Math.pow(time, 2.0)) + 
                   (initialVelocity * time);
    return displacement;
    }
   //solving for t
    public double kinematics14(double x, double a, double Vo){
        double t;
        double vF;
        vF = kinematics7(x, a, Vo);
        t = kinematics6(vF, a, Vo);
        return t;
    }
   //methods for second kinematic equation (v_f=v_o + at
   //solving for final velocity  
    public double kinematics3 (double initialVelocity, double time,
                               double acceleration) {
        double finalVelocity;
        finalVelocity = initialVelocity + (time * acceleration); 
        return finalVelocity;
    }
    //solving for initial velocity 
    public double kinematics4 (double finalVelocity, double time,
                               double acceleration){
        double initialVelocity;
        initialVelocity = finalVelocity - (time * acceleration);
        return initialVelocity;
    }
    //solving for acceleration
    public double kinematics5 (double finalVelocity, double time,
                               double initialVelocity) {
        double acceleration;
        acceleration = (finalVelocity - initialVelocity) / time;
        return acceleration;
    }
    //solving for time
    public double kinematics6 (double finalVelocity, double acceleration,
                               double initialVelocity) {
        double time;
        time = (finalVelocity - initialVelocity) / acceleration;
        return time;
    }
    //methods for third kinematic equation (v_f^2=v_o^2 + 2ax)
    //solving for final velocity
   public double kinematics7 (double displacement, double acceleration,
                               double initialVelocity) {
        double finalVelocity;
        finalVelocity = Math.sqrt(Math.pow(initialVelocity, 2.0) + 
                            (2 * acceleration * displacement));
        return finalVelocity;
    }
   //solving for initial velocity
   public double kinematics8 (double displacement, double acceleration,
                               double finalVelocity){
       double initialVelocity;
       initialVelocity = Math.sqrt(Math.pow(finalVelocity, 2.0) - 
                        (2 * acceleration * displacement));
       return initialVelocity;
   }
   // solving for displacement
   public double kinematics9 (double initialVelocity, double acceleration,
                               double finalVelocity){
       double displacement;
       displacement = (Math.pow(finalVelocity, 2.0) - Math.pow(initialVelocity, 2.0)) 
                       / (2 * acceleration);
       return displacement;
   }
   // solving for acceleration
   public double kinematics10 (double initialVelocity, double displacement,
                               double finalVelocity){
       double acceleration;
       acceleration = (Math.pow(finalVelocity, 2.0) - Math.pow(initialVelocity, 2.0)) 
                       / (2 * displacement);
       return acceleration;
   }
   //methods for fourth kinematic equation (x=v_ft-1/2at^2)
   //solving for displacement
   public double kinematics11 (double vf, double a, double t){
       double x;
       x = (vf * t) - (0.5 * a * t * t);
       return x;
   }
   //solving for v_f
   public double kinematics12 (double x, double a, double t){
       double vf;
       vf = (x + (0.5 * a * t * t)) / t;
       return vf;
   }
   //solving for a
   public double kinematics13 (double x, double t, double vf){
       double a;
       a = ((x - (vf * t)) * (-2.0)) / (t * t);
       return a;
   }
   //solving for t
   public double kinematics15 (double x, double vF, double a){
       double t;
       double vO;
       vO = kinematics8(x, a, vF);
       t = kinematics6(vF, a, vO);
       return t;
   } 
  /*
   * Method set to convert linear kinematics values into angular values
   */
  //Converting displacement to theta
   public double conversion1 (double radius, double displacement){
   double theta;
   theta = displacement / radius;
   return theta;
   }
  //Converting theta to displacement
   public double conversion2 (double radius, double theta){
   double displacement;
   displacement = radius * theta;
   return displacement;
   }
  //Converting velocity to sigma
   public double conversion3 (double radius, double velocity){
       double sigma;
       sigma = velocity / radius;
       return sigma;
   }
   //Converting sigma to velocity
   public double conversion4 (double sigma, double radius){
       double velocity;
       velocity = sigma * radius;
       return velocity;
   }
   //Converting alpha to acceleration
   public double conversion5 (double alpha, double radius){
       double acceleration;
       acceleration = alpha * radius;
       return acceleration;
   }
   //Converting acceleration to alpha
   public double conversion6 (double acceleration, double radius){
       double alpha;
       alpha = acceleration / radius; 
       return alpha;
   }
   /**
    * A set of methods to handle formulas involving Newton's Second Law
    * Can be used for both linear and angular applications, where
    * netForce = sumTorques
    * mass = moment of intertia (I)
    * acceleration = alpha
    */
   //Solving for netForce
   public double newton2nd1 (double acceleration, double mass){
       double netForce;
       netForce = mass * acceleration;
       return netForce;
   }
   //solving for mass or acceleration, depending on which parameter is entrered
   //can swith mass with acceleration to find the opposite one
   public double newton2n2 (double netForce, double mass){
       double acceleration;
       acceleration = netForce / mass;
       return acceleration;
   }
   //Methods to solve Knetic Energy equations ke=1/2mv^2 or ke=1/2Iw^2
   //Solve for ke
   public double kel(double mass, double v){
      double ke;
      ke = 0.5*mass*v*v;
        return ke;
   }
   //solve for m or I
   public double kelmass(double energy, double v){
       double mass;
       mass = 2*energy/v/v;
       return mass;
    }
   //solve for v or w
   public double kelv(double energy, double mass){
       double v;
       v = Math.sqrt(2*energy/mass);
       return v;
   }
   //Mehtods to solve equations of I
   //point mass
   public double MoI(double in[][], int numMasses){
       double moI;
       moI = 0.0;
       for(int i = 0; i<in.length; i++){
       moI=moI+(in[i][0]*in[i][1]);
       }
       return moI;
   }
   //Hoop around central axis
   public double moIHc(double mass, double radius){
   double moI=mass*radius*radius;
   return moI;
   }
   //Ring around central axis
   public double moIRi(double mass, double radius1, double radius2){
   double moI=.5*mass*(radius1*radius1+radius2*radius2);
   return moI;
   }
   //Cylinder around central axis
   public double moICc(double mass, double radius){
   double moI=.5*mass*radius*radius;
   return moI;
   }
   //Cylinder around central diameter
   public double moICd(double mass, double radius, double length){
   double moI=.25*mass*radius*radius+1/12*mass*length*length;
   return moI;
   }
   //Rod around central axis
   public double moIR(double mass, double length){
   double moI=1/12*mass*length*length;
   return moI;
   }
   //Solid Sphere around central axis
   public double moIS(double mass, double radius){
   double moI=2/5*mass*radius*radius;
   return moI;
   }
   //Hollow Sphere around central axis
   public double moISh(double mass, double radius){
   double moI=2/3*mass*radius*radius;
   return moI;
   }
   //Hoop around diameter
   public double moIHd(double mass, double radius){
   double moI=.5*mass*radius*radius;
   return moI;
   }
   //Slab around central axis
   public double moISl(double mass, double a, double b){
   double moI=mass/12*(a*a*b*b);
   return moI;
   }
   //Method to find additional I through the parallel axis theorum
   public double pAxis(double mass, double h){
   double mh2 = mass*h*h;
   return mh2;
   }
    public boolean inTest(String[] args) { 
    
    Scanner sc = new Scanner(System.in);
    String s = sc.next();
    try{
    Double.parseDouble(s);
    
    }catch(NumberFormatException e){
    return false;
    }
    return true;
    }
    
  }

