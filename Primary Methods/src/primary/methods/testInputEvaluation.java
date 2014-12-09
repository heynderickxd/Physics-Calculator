package primary.methods;
/**
 * @author heynderickxd
 */
import java.util.Scanner;
public class testInputEvaluation {
   public static void main(String[] args){
       //use a scanner to test if a user input is equal or not to a double
       Scanner input = new Scanner(System.in);
       String userType = input.next("0.987");
       if (input.hasNextDouble()){
           System.out.println("You have entered a correct input!");
       }
       else {
           System.out.println("You have not entered a correct input!");
       }
   } 
}
