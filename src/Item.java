/**
 * Represents a physical object and its attributes.
 *
 * @author Keeler Spear
 * @version %I%, %G%
 * @since 1.0
 */
public class Item {
    private static final int MIN_ROWS = 1;
    private static final int MIN_COLS = 1;
    private static final double tol = 0.00000000001;

    private double mass; //Mass (kg)
    private double area; //Cross-sectional area (m^2)
    private double drag; //Coefficient of drag

    public Item(double mass, double area, double drag) {
        this.mass = mass;
        this.area = area;
        this.drag = drag;
    }

    public double getMass() {
        return mass;
    }

    public double getArea() {
        return area;
    }

    public double getDrag() {
        return drag;
    }

    public void setMass(double mass) {
        this.mass = mass;
    }

    public void setArea(double area) {
        this.area = area;
    }

    public void setDrag(double drag) {
        this.drag = drag;
    }
}


