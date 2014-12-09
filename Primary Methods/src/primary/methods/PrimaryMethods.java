package primary.methods;
/**
 *Swag License 
 * @author heynderickxd
 * @author vannc
 * I used doubles all the way through. Must leave note reminding about sig figs 
 */
import java.util.Scanner;
public class PrimaryMethods {
   /**set of methods for dealing with kinematic equations
    * Note that the methods can be used for both linear and angular kinematics
    * problems. Simply substitute like this:
    * displacement = theta
    * initial velocity = sigma_o (w thing?)
    * final velocity = sigma_f 
    * acceleration = alfa

    */ 
   //methods for first kinematic equation (x=1/2*a*t^2 + v_ot)
   // solving for acceleration 
    public static double idk(String a, String b, String c, String d, String[] v){
        double r=0;
        for(int i=0; i<v.length; i++){
        if(d.equals(v[i])){
        r=r+i;
        }
        if(c.equals(v[i])){
        r=r+(10*i);
        }
        if(b.equals(v[i])){
        r=r+(100*i);
        }
        if(a.equals(v[i])){
        r=r+(1000*i);
        }
        }
        return r;
    }
    
    public static double kinematics1 (double displacement, double time, double initialVelocity) { // 431 341 4031 4301 3401 3041
    double acceleration = 0.0;
    acceleration = (displacement - (initialVelocity * time)) / Math.pow(time,2.0);
    return acceleration;
    }
   // solving for displacement
    public static double kinematics2 (double acceleration, double time, //1430 1340 4130 4310 3140 3410
                                double initialVelocity) {
    double displacement = 0.0;
    displacement = (0.5 * acceleration * Math.pow(time, 2.0)) + 
                   (initialVelocity * time);
    return displacement;
    }
       //solving for v_o
   public static double kinematics15 (double x, double a, double t){ //1403 1043 4013 4103 0143 0413
       double vo = 0.0;
       vo = (x - (0.5 * a * t * t)) / t;
       return vo;
   }
   //methods for second kinematic equation (v_f=v_o + at
   //solving for final velocity  
    public static double kinematics3 (double initialVelocity, double time, //3412 3142 4132 4312 1432 1342
                               double acceleration) {
        double finalVelocity = 0.0;
        finalVelocity = initialVelocity + (time * acceleration); 
        return finalVelocity;
    }
    //solving for initial velocity 
    public static double kinematics4 (double finalVelocity, double time, //2143 2413 4123 4213 1243 1423
                               double acceleration){
        double initialVelocity = 0.0;
        initialVelocity = finalVelocity - (time * acceleration);
        return initialVelocity;
    }
    //solving for acceleration
    public static double kinematics5 (double finalVelocity, double time, //3421 3241 4321 4231 2341 2431
                               double initialVelocity) {
        double acceleration = 0.0;
        acceleration = (finalVelocity - initialVelocity) / time;
        return acceleration;
    }
    //solving for time
    public static double kinematics6 (double finalVelocity, double acceleration, //3214 3124 1234 1324 2134 2314
                               double initialVelocity) {
        double time = 0.0;
        time = (finalVelocity - initialVelocity) / acceleration;
        return time;
    }
    //methods for third kinematic equation (v_f^2=v_o^2 + 2ax)
    //solving for final velocity
   public static double kinematics7 (double displacement, double acceleration, //1032 1302 3102 3012 0132 0312
                               double initialVelocity) {
        double finalVelocity = 0.0;
        finalVelocity = Math.sqrt(Math.pow(initialVelocity, 2.0) + 
                            (2 * acceleration * displacement));
        return finalVelocity;
    }
   //solving for initial velocity
   public static double kinematics8 (double displacement, double acceleration, //1023 1203 2013 2103 0123 0213
                               double finalVelocity){
       double initialVelocity = 0.0;
       initialVelocity = Math.sqrt(Math.pow(finalVelocity, 2.0) - 
                        (2 * acceleration * displacement));
       return initialVelocity;
   }
   // solving for displacement
   public static double kinematics9 (double initialVelocity, double acceleration, //3120 3210 2310 2130 1230 1320
                               double finalVelocity){
       double displacement = 0.0;
       displacement = (Math.pow(finalVelocity, 2.0) - Math.pow(initialVelocity, 2.0)) 
                       / (2 * acceleration);
       return displacement;
   }
   // solving for acceleration
   public static double kinematics10 (double initialVelocity, double displacement, //3021 3201 2031 2301 0231 0321
                               double finalVelocity){
       double acceleration = 0.0;
       acceleration = (Math.pow(finalVelocity, 2.0) - Math.pow(initialVelocity, 2.0)) 
                       / (2 * displacement);
       return acceleration;
   }
   //methods for fourth kinematic equation (x=v_ft-1/2at^2)
   //solving for displacement
   public static double kinematics11 (double vf, double a, double t){ //1420 1240 4210 4120 2410 2140
       double x = 0.0;
       x = (vf * t) - (0.5 * a * t * t);
       return x;
   }
   //solving for v_f
   public static double kinematics12 (double x, double a, double t){ //1402 1042 4012 4102 0142 0412
       double vf = 0.0;
       vf = (x + (0.5 * a * t * t)) / t;
       return vf;
   }
   //solving for a
   public static double kinematics13 (double x, double t, double vf){ //2401 2041 4201 4021 0241 0421
       double a = 0.0;
       a = ((x - (vf * t)) * (-2.0)) / (t * t);
       return a;
   }
  /**
   * Method set to convert linear kinematics values into angular values
   */
  //Converting displacement to theta
   public static double conversion1 (double radius, double displacement){
   double theta = 0.0;
   theta = displacement / radius;
   return theta;
   }
  //Converting theta to displacement
   public static double conversion2 (double radius, double theta){
   double displacement = 0.0;
   displacement = radius * theta;
   return displacement;
   }
  //Converting velocity to sigma
   public static double conversion3 (double radius, double velocity){
       double sigma = 0.0;
       sigma = velocity / radius;
       return sigma;
   }
   //Converting sigma to velocity
   public static double conversion4 (double sigma, double radius){
       double velocity = 0.0;
       velocity = sigma * radius;
       return velocity;
   }
   //Converting alpha to acceleration
   public static double conversion5 (double alpha, double radius){
       double acceleration = 0.0;
       acceleration = alpha * radius;
       return acceleration;
   }
   //Converting acceleration to alpha
   public static double conversion6 (double acceleration, double radius){
       double alpha = 0.0;
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
   public static double newton2nd1 (double acceleration, double mass){
       double netForce = 0.0;
       netForce = mass * acceleration;
       return netForce;
   }
   //solving for mass or acceleration, depending on which parameter is entrered
   //can swith mass with acceleration to find the opposite one
   public static double newton2n2 (double netForce, double mass){
       double acceleration = 0.0;
       acceleration = netForce / mass;
       return acceleration;
   }
   //Methods to solve Knetic Energy equations ke=1/2mv^2 or ke=1/2Iw^2
   //Solve for ke
   public static double kel(double mass, double v){
      double ke = 0.0;
      ke = 0.5*mass*v*v;
        return ke;
   }
   //solve for m or I
   public static double kelmass(double energy, double v){
       double mass = 0.0;
       mass = 2*energy/v/v;
       return mass;
    }
   //solve for v or w
   public static double kelv(double energy, double mass){
       double v = 0.0;
       v = Math.sqrt(2*energy/mass);
       return v;
   }
   //Mehtods to solve equations of I
   //point mass
   public static double MoI(double in[][], int numMasses){
       double moI = 0.0;
       for(int i = 0; i<in.length; i++){
       moI=moI+(in[i][0]*in[i][1]);
       }
       return moI;
   }
   //Hoop around central axis
   public static double moIHc(double mass, double radius){
   double moI=mass*radius*radius;
   return moI;
   }
   //Ring around central axis
   public static double moIRi(double mass, double radius1, double radius2){
   double moI=.5*mass*(radius1*radius1+radius2*radius2);
   return moI;
   }
   //Cylinder around central axis
   public static double moICc(double mass, double radius){
   double moI=.5*mass*radius*radius;
   return moI;
   }
   //Cylinder around central diameter
   public static double moICd(double mass, double radius, double length){
   double moI=.25*mass*radius*radius+1/12*mass*length*length;
   return moI;
   }
   //Rod around central axis
   public static double moIR(double mass, double length){
   double moI=1/12*mass*length*length;
   return moI;
   }
   //Solid Sphere around central axis
   public static double moIS(double mass, double radius){
   double moI=2/5*mass*radius*radius;
   return moI;
   }
   //Hollow Sphere around central axis
   public static double moISh(double mass, double radius){
   double moI=2/3*mass*radius*radius;
   return moI;
   }
   //Hoop around diameter
   public static double moIHd(double mass, double radius){
   double moI=.5*mass*radius*radius;
   return moI;
   }
   //Slab around central axis
   public static double moISl(double mass, double a, double b){
   double moI=mass/12*(a*a*b*b);
   return moI;
   }
   //Method to find additional I through the parallel axis theorum
   public static double pAxis(double mass, double h){
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

