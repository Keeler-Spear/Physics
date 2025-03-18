import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

public class Orbit {

    final static double GM = 4 * Math.PI * Math.PI; //Grav. const. * Mass of Sun (au^3/yr^2)
    final static double MASS = 1; //Mass of the commit
    final static double ADAPTRK = 1e-3;

    //x = [t, rx, ry, vx, vy, ax, ay]
    static final NFunction<Double> rx = (x) -> x[3];
    static final NFunction<Double> ry = (x) -> x[4];
    static final NFunction<Double> vx = (x) -> -GM * x[1] / Math.pow(Math.sqrt(Math.pow(x[1], 2) + Math.pow(x[2], 2)), 3);
    static final NFunction<Double> vy = (x) -> -GM * x[2] / Math.pow(Math.sqrt(Math.pow(x[1], 2) + Math.pow(x[2], 2)), 3);

    static NFunction<Double>[] system = new NFunction[]{rx, ry, vx, vy};

    /**
     * Calculates the orbit of a commit.
     *
     * @param r0 The initial radial distance, in AU. 2x1
     * @param v0 The initial tangential velocity of the projectile, in AU/yr. 2x1
     * @param numStep The number of steps the program should take
     * @param h The time step for the program to use, in yr.
     * @param method The numerical method to be used. 1 = Euler, 2 = Euler-Cromer, 3 = RK4, 4 = Adaptive RK4, 5 = AB4
     * @return An array of matrices containing:
     *         <ul>
     *             <li> A matrix containing the angle values of the motion</li>
     *             <li> A matrix containing the radius value of the motion</li>
     *             <li> A matrix containing the time values corresponding to the angles and radius of motion</li>
     *             <li> A Matrix containing the kinetic energy of the comet at each time step</li>
     *             <li> A Matrix containing the potential energy of the comet at each time step</li>
     *         </ul>
     * @throws IllegalArgumentException If the method selected is invalid.
     */
    public static Matrix[] calculate(double r0, double v0, int numStep, double h, int method) {
        if (method > 5) {
            throw new IllegalArgumentException("Invalid Method!");
        }

        Matrix r = new Matrix(new double[] {r0, 0.0});
        Matrix v = new Matrix(new double[] {0.0, v0});
        Matrix a;

        Matrix rPlot = new Matrix (2, 1);
        Matrix thPlot = new Matrix (2, 1);
        Matrix tPlot = new Matrix (2, 1);
        Matrix realTPlot = new Matrix (2, 1);
        Matrix kinetic = new Matrix (2, 1);
        Matrix potential = new Matrix(numStep, 1);
        Matrix sol = new Matrix(2, 1);
        double t = 0.0;
        double rNorm;
        double[] initialConditions = {r0, 0.0, 0.0, v0};

        if (method != 1 && method != 2) { //Solve the ODE using an external method
            if (method == 3) { //RK4
                sol = ODE.rk4System(system, t, initialConditions, t + h * numStep, h); //[rx, ry, vx, vy, ax, ay]
            }
            else if (method == 4) { //Adaptive RK4
                sol = ODE.adaptiveRK4System(system, t, initialConditions, t + h * numStep, h, 0.9, 4.0);
                realTPlot = LinearAlgebra.vectorFromColumn(sol, sol.getCols());
                sol.removeCol(sol.getCols()); //Removing the time values
            }
            else { //AB4
                throw new IllegalArgumentException("AB4 is not yeet supported!");
            }
        }

        if (method == 1 || method == 2 || method == 3) {

            rPlot = new Matrix(numStep, 1);
            thPlot = new Matrix(numStep, 1);
            tPlot = new Matrix(numStep, 1);
            kinetic = new Matrix(numStep, 1);
            potential = new Matrix(numStep, 1);

            for (int i = 1; i <= numStep; i++) {
                //Recording the current values
                rNorm = LinearAlgebra.l2Norm(r);
                rPlot.setValue(i, 1, rNorm);
                thPlot.setValue(i, 1, Math.atan2(r.getValue(2, 1), r.getValue(1, 1)));
                tPlot.setValue(i, 1, t);
                kinetic.setValue(i, 1, 0.5 * MASS * Math.pow(LinearAlgebra.l2Norm(v), 2));
                potential.setValue(i, 1, -GM * MASS / rNorm);

                if (method == 1) { //Euler
                    a = LinearAlgebra.scaleMatrix(r, -GM / Math.pow(rNorm, 3));
                    r = LinearAlgebra.addMatrices(r, v, h);
                    v = LinearAlgebra.addMatrices(v, a, h);
                } else if (method == 2) { //Euler-Cromer
                    a = LinearAlgebra.scaleMatrix(r, -GM / Math.pow(rNorm, 3));
                    v = LinearAlgebra.addMatrices(v, a, h);
                    r = LinearAlgebra.addMatrices(r, v, h);
                } else { //Other methods
                    r.setValue(1, 1, sol.getValue(i, 1));
                    r.setValue(2, 1, sol.getValue(i, 2));
                    v.setValue(1, 1, sol.getValue(i, 3));
                    v.setValue(2, 1, sol.getValue(i, 4));
                }

                t += h;
            }
        }

        if (method == 4) {
            tPlot = realTPlot;
            rPlot = new Matrix(sol.getRows(), 1);
            thPlot = new Matrix(sol.getRows(), 1);
            kinetic = new Matrix(sol.getRows(), 1);
            potential = new Matrix(sol.getRows(), 1);



            for (int i = 1; i <= sol.getRows(); i++) {
                rNorm = LinearAlgebra.l2Norm(r);
                rPlot.setValue(i, 1, rNorm);
                thPlot.setValue(i, 1, Math.atan2(r.getValue(2, 1), r.getValue(1, 1)));
                kinetic.setValue(i, 1, 0.5 * MASS * Math.pow(LinearAlgebra.l2Norm(v), 2));
                potential.setValue(i, 1, -GM * MASS / rNorm);

                r.setValue(1, 1, sol.getValue(i, 1));
                r.setValue(2, 1, sol.getValue(i, 2));
                v.setValue(1, 1, sol.getValue(i, 3));
                v.setValue(2, 1, sol.getValue(i, 4));
            }
        }

        return new Matrix[] {thPlot, rPlot, tPlot, kinetic, potential};
    }


    public static void plotTrajectory(Matrix theta, Matrix r) {
        try {
            File pythonScript = File.createTempFile("trajectory_plot", ".py");
            pythonScript.deleteOnExit();

            try (PrintWriter out = new PrintWriter(new FileWriter(pythonScript))) {
                //Imports
                out.println("import matplotlib.pyplot as plt");
                out.println("import numpy as np");
                //Initializing data
                out.println("thplot = np.array(" + theta.npString() + ")");
                out.println("rplot = np.array(" + r.npString() + ")");
                //Creating the plot
                out.println("ax = plt.subplot(111, projection='polar')");
                out.println("ax.plot(thplot,rplot,'+')");
                out.println("ax.set_title('Distance (AU)')");
                out.println("ax.grid(True)");
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

    public static void plotEnergy(Matrix tPlot, Matrix kinetic, Matrix potential) {
        Matrix totalEnergy = LinearAlgebra.addMatrices(kinetic, potential, 1);

        try {
            File pythonScript = File.createTempFile("energy_plot", ".py");
            pythonScript.deleteOnExit();

            try (PrintWriter out = new PrintWriter(new FileWriter(pythonScript))) {
                //Imports
                out.println("import matplotlib.pyplot as plt");
                out.println("import numpy as np");
                //Initializing data
                out.println("tPlot = np.array(" + tPlot.npString() + ")");
                out.println("kinetic = np.array(" + kinetic.npString() + ")");
                out.println("potential = np.array(" + potential.npString() + ")");
                out.println("totalE = np.array(" + totalEnergy.npString() + ")");
                //Creating the plot
                out.println("plt.plot(tPlot,kinetic,'-.',tPlot,potential,'--',tPlot,totalE,'-')");
                out.println("plt.legend(['Kinetic','Potential','Total']);");
                out.println("plt.xlabel('Time (yr)')");
                out.println("plt.ylabel(r'Energy ($M AU^2/yr^2$)')");
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

    public static void plotStepVsV0(Matrix hs, Matrix vs) {

        try {
            File pythonScript = File.createTempFile("stepvv_plot", ".py");
            pythonScript.deleteOnExit();

            try (PrintWriter out = new PrintWriter(new FileWriter(pythonScript))) {
                //Imports
                out.println("import matplotlib.pyplot as plt");
                out.println("import numpy as np");
                //Initializing data
                out.println("hplot = np.array(" + hs.npString() + ")");
                out.println("vplot = np.array(" + vs.npString() + ")");
                //Creating the plot
                out.println("plt.loglog(vplot, hplot, marker='o', linestyle='-', color='b', label='Given Data')");
                out.println("plt.xlabel('Initial Velocity AU/yr')");
                out.println("plt.ylabel('Step Size')");
                out.println("plt.title('Step Size VS. Initial Velocity')");
                out.println("plt.legend()");
                out.println("plt.grid(True, which='both', linestyle='--', linewidth=0.5)");
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

}
