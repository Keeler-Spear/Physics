import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;

public class Simulation {
    final static int MAX_ITERATIONS = 1000;
    final static Matrix g = new Matrix(new double[] {0, -9.806}); //Gravitational constant (m/s^2)
    final static Matrix noWind = new Matrix(new double[] {0.0, 0.0});
    final static Matrix vDrag = new Matrix(new double[] {0, 11.176, 22.352, 33.528, 44.704, 55.88});
    final static Matrix cDrag = new Matrix(new double[] {0.5, 0.5, 0.5, 0.4, 0.28, 0.23});

    /**
     * Calculates the motion of a projectile.
     *
     * @param height The height the projectile will be launched from, in meters.
     * @param speed The initial speed of the projectile, in meters per second.
     * @param angle The angle the projectile in launched at, in degrees.
     * @param timeStep The time step for the program to use.
     * @param method The numerical method to be used. 1 = Euler, 2 = Euler-Cromer, 3 = Midpoint.
     * @param air If air resistance should be factored into the motion
     * @param wind The speed of the wind.
     * @param constantDrag If the drag coefficient of the object is constant.
     * @return A matrix containing the x and y positions of the projectile at each time step.
     */
    public static Matrix projectileMotion(Item object, double height, double speed, double angle, double timeStep, int method, boolean air, Matrix wind, boolean constantDrag) {

        //Initializing Variables
        double t = 0; //Current time (seconds)
        ArrayList<Double> time = new ArrayList<>(); //Array list to store time values
        ArrayList<Double> x = new ArrayList<>(); //Array list to store x values
        ArrayList<Double> y = new ArrayList<>(); //Array list to store y values

        double mass = object.getMass(); //Mass (km)
        double area = object.getArea(); //Cross-sectional area (m^2)
        double coD = object.getDrag(); //Coefficient of drag
        double density = 1.2; //Density of air (kg/m^3)
        double airRes;
        if (air) {
            airRes = -0.5 * coD * density * area / mass; //Air Resistance
        }
        else {
            airRes = 0.0; //If we are not using air resistance, set air resistance to zero
        }

        Matrix changingCOD = null;

        if (!constantDrag) {
            changingCOD = Interpolation.polynomial(vDrag, cDrag, LinearAlgebra.linSpace(0, 100, timeStep)); //If drag is not constant, interpolate drag values
        }

        Matrix r = new Matrix(new double[] {0.0, height}); //Position vector
        Matrix v = new Matrix(new double[] {speed * Math.cos(angle * Math.PI / 180), speed * Math.sin(angle * Math.PI / 180)}); //Velocity vector
        Matrix vRel;
        Matrix a; //Acceleration vector
        Matrix windAcc = LinearAlgebra.scaleMatrix(wind, -airRes * LinearAlgebra.l2Norm(wind)); //Wind = airRes * |wind| * wind
        int j = 0;

        //Calculations
        if (method == 1) { //Euler
            while ((r.getValue(2, 1) > 0 || j == 0) && j < MAX_ITERATIONS) {
                x.add(r.getValue(1, 1));
                y.add(r.getValue(2, 1));
                time.add(t);
                t += timeStep;

                vRel = LinearAlgebra.addMatrices(v, wind, 1);

                r = LinearAlgebra.addMatrices(r, vRel, timeStep);

                if (!constantDrag) {
                    airRes = -0.5 * changingCOD.getValue(j+1, 1) * density * area / mass;
                }

                a = LinearAlgebra.scaleMatrix(vRel, airRes *  LinearAlgebra.l2Norm(vRel));
                a = LinearAlgebra.addMatrices(a, g, 1);

                v = LinearAlgebra.addMatrices(v, a, timeStep);
                j++;
            }
        }

        if (method == 2) { //Euler-Cromer
            while ((r.getValue(2, 1) > 0 || j == 0) && j < MAX_ITERATIONS) {
                x.add(r.getValue(1, 1));
                y.add(r.getValue(2, 1));
                time.add(t);
                t += timeStep;

                if (!constantDrag) {
                    airRes = -0.5 * changingCOD.getValue(j+1, 1) * density * area / mass;
                }

                vRel = LinearAlgebra.addMatrices(v, wind, 1);
                a = LinearAlgebra.scaleMatrix(vRel, airRes *  LinearAlgebra.l2Norm(vRel));
                a = LinearAlgebra.addMatrices(a, g, 1);

                v = LinearAlgebra.addMatrices(v, a, timeStep);
                vRel = LinearAlgebra.addMatrices(v, wind, 1);

                r = LinearAlgebra.addMatrices(r, vRel, timeStep);

                j++;
            }
        }

        else { //Midpoint
            while ((r.getValue(2, 1) > 0 || j == 0) && j < MAX_ITERATIONS) {
                x.add(r.getValue(1, 1));
                y.add(r.getValue(2, 1));
                time.add(t);
                t += timeStep;


                if (!constantDrag) {
                    airRes = -0.5 * changingCOD.getValue(j+1, 1) * density * area / mass;
                }

                vRel = LinearAlgebra.addMatrices(v, wind, 1);
                a = LinearAlgebra.scaleMatrix(vRel, airRes *  LinearAlgebra.l2Norm(vRel));
                a = LinearAlgebra.addMatrices(a, g, 1);

                r = LinearAlgebra.addMatrices(r, LinearAlgebra.addMatrices(vRel, a, timeStep / 2), timeStep);

                v = LinearAlgebra.addMatrices(v, a, timeStep);

                j++;
            }
        }

        x.add(r.getValue(1, 1));
        y.add(r.getValue(2, 1));
        time.add(t);

        Matrix vals = new Matrix (x.size(), 3);
        for (int i = 0; i < x.size(); i++) {
            vals.setValue(i + 1, 1, time.get(i));
            vals.setValue(i + 1, 2, x.get(i));
            vals.setValue(i + 1, 3, y.get(i));
        }


        return vals;
    }

    //a = airRes * |vRel|vRel + g
    private static Matrix calcWindAcc(Matrix v, Matrix wind, double airRes) {
        Matrix vRel = LinearAlgebra.addMatrices(v, wind, -1);
        Matrix a = LinearAlgebra.scaleMatrix(vRel, airRes * LinearAlgebra.l2Norm(vRel));
        a = LinearAlgebra.addMatrices(a, g, 1);
        return a;
    }

    public static void displayMotion(Matrix theory, Matrix airRes) {
        Matrix xNoAir = LinearAlgebra.vectorFromColumn(theory, 2);
        Matrix yNoAir = LinearAlgebra.vectorFromColumn(theory, 3);

        Matrix x = LinearAlgebra.vectorFromColumn(airRes,2);
        Matrix y = LinearAlgebra.vectorFromColumn(airRes,3);


        int lastStep = airRes.getRows();

        try {
            File pythonScript = File.createTempFile("scatter_plot", ".py");
            pythonScript.deleteOnExit();

            try (PrintWriter out = new PrintWriter(new FileWriter(pythonScript))) {
                //Imports
                out.println("import matplotlib.pyplot as plt");
                out.println("import numpy as np");
                //Initializing data
                out.println("x = np.array(" + x.npString() + ")");
                out.println("y = np.array(" + y.npString() + ")");
                out.println("xNoAir = np.array(" + xNoAir.npString() + ")");
                out.println("yNoAir = np.array(" + yNoAir.npString() + ")");
                out.println("xground = np.array([0., " + xNoAir.getValue(lastStep, 1) + "])");
                out.println("yground = np.array([0., 0.])");
                out.println("plt.plot(x, y, '+', xNoAir[0:" + lastStep + "], yNoAir[0:" + lastStep + "], '-', xground,yground,'r-')");
                //Printing the scatter plot
                out.println("plt.legend(['Euler method', 'Theory (No air)']);");
                out.println("plt.xlabel('Range (m)')");
                out.println("plt.ylabel('Height (m)')");
                out.println("plt.title('Projectile motion')");
                out.println("plt.show()");
            }

            ProcessBuilder pb = new ProcessBuilder("python", pythonScript.getAbsolutePath());
            pb.inheritIO();
            Process p = pb.start();
            p.waitFor();
        } catch (IOException | InterruptedException e) {
            e.printStackTrace();
        }
    }

    public static double flightStats(Matrix r, boolean printComparison, boolean printStats) {
        Matrix x = LinearAlgebra.vectorFromColumn(r, 2);
        Matrix y = LinearAlgebra.vectorFromColumn(r, 3);

        Matrix lastThreeTVals = new Matrix(new double[] {r.getValue(r.getRows() - 2, 1), r.getValue(r.getRows() - 1, 1), r.getValue(r.getRows(), 1)});
        Matrix lastThreeXVals = new Matrix(new double[] {x.getValue(x.getRows() - 2, 1), x.getValue(x.getRows() - 1, 1), x.getValue(x.getRows(), 1)});
        Matrix lastThreeYVals = new Matrix(new double[] {y.getValue(y.getRows() - 2, 1), y.getValue(y.getRows() - 1, 1), y.getValue(y.getRows(), 1)});

        Matrix xIntercept = Interpolation.polynomial(lastThreeYVals, lastThreeXVals, new Matrix(new double[] {0}));

        double maxRange = x.getValue(x.getRows(), 1);
        double maxRangeCorrected = xIntercept.getValue(1, 1);
        Matrix tIntercept = Interpolation.polynomial(lastThreeXVals, lastThreeTVals, new Matrix(new double[] {maxRangeCorrected}));
        double time = r.getValue(r.getRows(), 1);
        double timeCorrected = tIntercept.getValue(1, 1);

        if (printComparison) {
            System.out.print("Uncorrected maximum range is " + String.format("%.3f", maxRange) + " meters");
            System.out.print("\tCorrected maximum range is " + String.format("%.3f", maxRangeCorrected) + " meters");
            System.out.println("\tAbsolute difference: " + String.format("%.3f", Error.absolute(maxRange, maxRangeCorrected)) + " meters");
            System.out.print("Uncorrected time of flight is " + String.format("%.3f", time) + " seconds");
            System.out.print("\tCorrected time of flight is " + String.format("%.3f", timeCorrected) + " seconds");
            System.out.println("\tAbsolute difference: " + String.format("%.3f", Error.absolute(time, timeCorrected)) + " seconds");
        }

        else if (printStats) {
            System.out.println("Maximum range: " + String.format("%.3f", maxRangeCorrected) + " meters");
            System.out.println("Flight time: " + String.format("%.3f", timeCorrected) + " seconds");
        }
        return maxRangeCorrected;
    }

    public static void plotRangeVSAngle (Item object, double height, double speed, double timeStep, int method, boolean air, Matrix wind, boolean constantDrag, double a, double b) {
        Matrix angles = LinearAlgebra.linSpace(a, b, 0.1);
        Matrix ranges = new Matrix(angles.getRows(), 1);
        double range;
        for (int i = 1; i <= ranges.getRows(); i++) {
            range = flightStats(projectileMotion(object, height, speed, angles.getValue(i, 1), timeStep, method, air, wind, constantDrag), false, false);
            ranges.setValue(i, 1, range);
        }

        int index = 1;
        double max = 0.0;
        for (int i = 1; i <= ranges.getRows(); i++) {
            if (max < ranges.getValue(i, 1)) {
                max = ranges.getValue(i, 1);
                index = i;
            }
        }
        System.out.println("Best Angle: " + angles.getValue(index, 1));

        PyChart.scatter(angles, ranges, "Projectile Range", "Angle (degrees)", "Range (meters)", "Range Vs. Angle");
    }

    public static void plotRangeVSWind (Item object, double height, double speed, double angle, double timeStep, int method, boolean constantDrag, double a, double b) {
        Matrix winds = LinearAlgebra.linSpace(a, b, 5);
        Matrix ranges = new Matrix(winds.getRows(), 1);
        double range;
        Matrix wind;
        for (int i = 1; i <= ranges.getRows(); i++) {
            wind = new Matrix(new double[] {winds.getValue(i, 1), 0.0});
            range = flightStats(projectileMotion(object, height, speed, angle, timeStep, method, true, wind, constantDrag), false, false);
            ranges.setValue(i, 1, range);
        }

        PyChart.fnc(winds, ranges, "Projectile Range", "Wind Velocity (m/s)", "Range (meters)", "Range Vs. Wind Velocity");
    }

    public static void plotTwoRangeVSAngle (Item object, double height, double speed, double timeStep, int method, boolean air, Matrix wind, double a, double b) {
        Matrix angles = LinearAlgebra.linSpace(a, b, 0.1);
        Matrix conRanges = new Matrix(angles.getRows(), 1);
        double range;
        for (int i = 1; i <= conRanges.getRows(); i++) {
            range = flightStats(projectileMotion(object, height, speed, angles.getValue(i, 1), timeStep, method, air, wind, true), false, false);
            conRanges.setValue(i, 1, range);
        }

        Matrix delRanges = new Matrix(angles.getRows(), 1);
        for (int i = 1; i <= conRanges.getRows(); i++) {
            range = flightStats(projectileMotion(object, height, speed, angles.getValue(i, 1), timeStep, method, air, wind, false), false, false);
            delRanges.setValue(i, 1, range);
        }

        PyChart.twoFnc(angles, conRanges, "Constant Drag", delRanges, "Changing Drag", "Angle (degrees)", "Range (meters)", "Range Vs. Angle");
    }

}
