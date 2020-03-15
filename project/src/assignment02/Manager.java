package assignment02;


import ij.IJ;
import ij.ImagePlus;
import ij.plugin.filter.PlugInFilter;
import ij.process.ImageProcessor;
import imagingbook.pub.geometry.basic.Point;
import imagingbook.pub.geometry.delaunay.Triangle;
import imagingbook.pub.geometry.delaunay.chew.TriangulationChew;
import imagingbook.pub.geometry.*;

//import java.awt.*;
//import java.awt.geom.Point2D;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

import imagingbook.pub.geometry.*;
import org.apache.commons.math3.linear.DecompositionSolver;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.linear.SingularValueDecomposition;

class Manager implements PlugInFilter {
    private TriangulationChew triangulator_X;
    private TriangulationChew triangulator_X_2;


    public int setup(String arg, ImagePlus im) {
        return PlugInFilter.DOES_8G;
    }

    public void run(ImageProcessor ip) {
        final int w = ip.getWidth();
        final int h = ip.getHeight();

        // Collect all image points with pixel values greater than zero:
        List<Point> points_X = new ArrayList<Point>();

        for (int v = 0; v < h; v++) {
            for (int u = 0; u < w; u++) {
                int p = ip.getPixel(u, v);
                if (p > 0) {
                    points_X.add(Point.create(u, v));
                }
            }
        }
        triangulator_X = new TriangulationChew(points_X);
        List<Triangle> triangles_X = triangulator_X.getTriangles();
        triangulator_X_2 = new TriangulationChew(points_X_2);


        for (int i = 0; i < triangles_X.size(); i++) {
            Point[] points_from_X = triangles_X.get(i).getPoints();
            Point[] points_from_X_prim = triangles_X.get(i).getPoints();
            double x_x0 = points_from_X[0].getX();
            double x_y0 = points_from_X[0].getY();
            double x_x1 = points_from_X[0].getX();
            double x_y1 = points_from_X[0].getY();
            double x_x2 = points_from_X[0].getX();
            double x_y2 = points_from_X[0].getY();
            double x2_x0 = points_from_X_prim[0].getX();
            double x2_y0 = points_from_X_prim[0].getY();
            double x2_x1 = points_from_X_prim[0].getX();
            double x2_y1 = points_from_X_prim[0].getY();
            double x2_x2 = points_from_X_prim[0].getX();
            double x2_y2 = points_from_X_prim[0].getY();
            RealMatrix A = MatrixUtils.createRealMatrix(new double[][]
                    {{x_x0,x_y0,1,0,0,0},
                    {0,0,0,x_x0,x_y0,1},
                    {x_x1,x_y1,1,0,0,0},
                    {0,0,0,x_x1,x_y1,1},
                    {x_x2,x_y2,1,0,0,0},
                    {0,0,0,x_x2,x_y2,1},
                    });
            RealVector b = MatrixUtils.createRealVector(new double[]{ x2_x0, x2_y0, x2_x1, x2_y1, x2_x2, x2_y2);
            DecompositionSolver s = new SingularValueDecomposition(A).getSolver();
            RealVector affine_transformation = s.solve(b);




        }
        }


    }