import java.util.ArrayList;

public class Simulation {
    final static double TOL = 0.000001;
    final static int MAX_ITERATIONS = 1000;

    //Currently models a baseball. Should add an input for an acceleration vector and an input for an "object"
    /**
     * Calculates the motion of a projectile.
     *
     * @param height The height the projectile will be launched from, in meters.
     * @param speed The initial speed of the projectile, in meters per second.
     * @param angle The angle the projectile in launched at, in degrees.
     * @param timeStep The time step for the program to use.
     * @return TODO
     * @throws IllegalArgumentException TODO
     */
    public static void projectileMotion(double height, double speed, double angle, double timeStep) {
        double t = 0;
        //n x 3 Matrix with [time, xVal, yVal]
        ArrayList<Double> x = new ArrayList<>();
        ArrayList<Double> y = new ArrayList<>();

        double coD = 0.35; //Coefficient of drag
        double area = 4.3e-3; //Cross-sectional area (m^2)
        double g = 9.81; //Gravitational constant (m/s^2)
        double mass = 0.145; //Mass (km)
        double density = 1.2; //Density of air (kg/m^3)
        double airRes = -0.5 * coD * density * area / mass; //Air Resistance

        Matrix r0 = new Matrix(new double[] {0.0, height});
        Matrix v0 = new Matrix(new double[] {speed * Math.cos(angle * Math.PI / 180), speed * Math.sin(angle * Math.PI / 180)});
        Matrix r = new Matrix(new double[] {0.0, height});
        Matrix v = new Matrix(new double[] {speed * Math.cos(angle * Math.PI / 180), speed * Math.sin(angle * Math.PI / 180)});
        Matrix a = new Matrix(new double[] {airRes, g});

        System.out.println(v0);
//
//        while (r[1] > 0) {
//            x.add(r[0]);
//            y.add(r[1]);
//            t += timeStep;
//
//            r[0] += timeStep * v[0];
//            r[1] += timeStep * v[1];
//
//            v[0] +=
//        }

    }

}
