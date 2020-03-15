package assignment02;


import ij.ImagePlus;
import ij.plugin.filter.PlugInFilter;
import ij.process.ImageProcessor;
import org.apache.commons.math3.linear.DecompositionSolver;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.linear.SingularValueDecomposition;


public class AffineTransformator
{
    private RealMatrix transformations;
    private RealMatrix transformed_points;

    public AffineTransformator(RealMatrix transformations)
    {
        this.transformations=transformations;
    }


    public RealMatrix getTransformed_points() {
        return transformed_points;
    }

    public void transform(RealMatrix points)
    {
        this.transformed_points=this.transformations.multiply(points);
    }
    public double calculate_error(RealMatrix points)
    {
        RealMatrix substracted=this.transformed_points.subtract(points);
        RealMatrix substracted_power=substracted.transpose().multiply(substracted);
        double error=0;
        for(int i=0;i<substracted_power.getColumnDimension();i++)
        {
            error+=Math.sqrt(substracted_power.getColumn(i)[0]+substracted_power.getColumn(i)[1]);
        }
        return error;
    }


}