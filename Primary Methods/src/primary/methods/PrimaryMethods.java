package primary.methods;
/**
 *
 * @author heynderickxd
 * I used doubles all the way through. Must leave note reminding about sig figs 
 */
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
    public double kinematics1 (double displacement, double time, 
                            double initialVelocity) {
    double acceleration = 0.0;
    acceleration = (displacement - (initialVelocity * time)) / Math.pow(time,2.0);
    return acceleration;
    }
   // solving for displacement
    public double kinematics2 (double acceleration, double time, 
                                double initialVelocity) {
    double displacement = 0.0;
    displacement = (0.5 * acceleration * Math.pow(time, 2.0)) + 
                   (initialVelocity * time);
    return displacement;
    }
   //methods for second kinematic equation (v_f=v_o + at
   //solving for final velocity  
    public double kinematics3 (double initialVelocity, double time,
                               double acceleration) {
        double finalVelocity = 0.0;
        finalVelocity = initialVelocity + (time * acceleration); 
        return finalVelocity;
    }
    //solving for initial velocity 
    public double kinematics4 (double finalVelocity, double time,
                               double acceleration){
        double initialVelocity = 0.0;
        initialVelocity = finalVelocity - (time * acceleration);
        return initialVelocity;
    }
    //solving for acceleration
    public double kinematics5 (double finalVelocity, double time,
                               double initialVelocity) {
        double acceleration = 0.0;
        acceleration = (finalVelocity - initialVelocity) / time;
        return acceleration;
    }
    //solving for time
    public double kinematics6 (double finalVelocity, double acceleration,
                               double initialVelocity) {
        double time = 0.0;
        time = (finalVelocity - initialVelocity) / acceleration;
        return time;
    }
    //methods for third kinematic equation (v_f^2=v_o^2 + 2ax)
    //solving for final velocity
   public double kinematics7 (double displacement, double acceleration,
                               double initialVelocity) {
        double finalVelocity = 0.0;
        finalVelocity = Math.sqrt(Math.pow(initialVelocity, 2.0) + 
                            (2 * acceleration * displacement));
        return finalVelocity;
    }
   //solving for initial velocity
   public double kinematics8 (double displacement, double acceleration,
                               double finalVelocity){
       double initialVelocity = 0.0;
       initialVelocity = Math.sqrt(Math.pow(finalVelocity, 2.0) - 
                        (2 * acceleration * displacement));
       return initialVelocity;
   }
   // solving for displacement
   public double kinematics9 (double initialVelocity, double acceleration,
                               double finalVelocity){
       double displacement = 0.0;
       displacement = (Math.pow(finalVelocity, 2.0) - Math.pow(initialVelocity, 2.0)) 
                       / (2 * acceleration);
       return displacement;
   }
   // solving for acceleration
   public double kinematics10 (double initialVelocity, double displacement,
                               double finalVelocity){
       double acceleration = 0.0;
       acceleration = (Math.pow(finalVelocity, 2.0) - Math.pow(initialVelocity, 2.0)) 
                       / (2 * displacement);
       return acceleration;
   }
  /**
   * Method set to convert linear kinematics values into angular values
   */
  //Converting displacement to theta
   public double conversion1 (double radius, double displacement){
   double theta = 0.0;
   theta = displacement / radius;
   return theta;
   }
  //Converting theta to displacement
   public double conversion2 (double radius, double theta){
   double displacement = 0.0;
   displacement = radius * theta;
   return displacement;
   }
  //Converting velocity to sigma
   public double conversion3 (double radius, double velocity){
       double sigma = 0.0;
       sigma = velocity / radius;
       return sigma;
   }
   //Converting sigma to velocity
   public double conversion4 (double sigma, double radius){
       double velocity = 0.0;
       velocity = sigma * radius;
       return velocity;
   }
   //Converting alpha to acceleration
   public double conversion5 (double alpha, double radius){
       double acceleration = 0.0;
       acceleration = alpha * radius;
       return acceleration;
   }
   //Converting acceleration to alpha
   public double conversion6 (double acceleration, double radius){
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
   public double newton2nd1 (double acceleration, double mass){
       double netForce = 0.0;
       netForce = mass * acceleration;
       return netForce;
   }
   //solving for mass or acceleration, depending on which parameter is entrered
   //can swith mass with acceleration to find the opposite one
   public double newton2n2 (double netForce, double mass){
       double acceleration = 0.0;
       acceleration = netForce / mass;
       return acceleration;
   }
}
