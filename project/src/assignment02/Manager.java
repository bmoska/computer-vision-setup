package assignment02;


import ij.IJ;
import ij.ImagePlus;
import ij.plugin.filter.PlugInFilter;
import ij.process.ByteProcessor;
import ij.process.ImageProcessor;
import ij.ImageStack;
import imagingbook.pub.geometry.basic.Point;
import imagingbook.pub.geometry.delaunay.Triangle;
import imagingbook.pub.geometry.delaunay.chew.TriangulationChew;

//import java.awt.*;
//import java.awt.geom.Point2D;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

import org.apache.commons.math3.linear.DecompositionSolver;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.linear.SingularValueDecomposition;

class Manager implements PlugInFilter {
    private TriangulationChew triangulator_X;
    private TriangulationChew triangulator_X_2;
    ImagePlus im = null;


    public int setup(String arg, ImagePlus im) {

        this.im = im;
        return DOES_8G + STACK_REQUIRED;
    }
    public void run(ImageProcessor ip) {

        ImageStack stack = im.getStack();

        // Collect all image points with pixel values greater than zero:
        //List<Point> points_X = new ArrayList<Point>();
        List<List<Point>> points = new ArrayList<List<Point>>();
        for (int k = 1; k <= 2; k++) {    // NOTE: k starts from 1!!
            ImageProcessor ipk = stack.getProcessor(k);
            List<Point> points_of_single_im = new ArrayList<Point>();
            final int w = ipk.getWidth();
            final int h = ipk.getHeight();
            for (int v = 0; v < h; v++) {
                for (int u = 0; u < w; u++) {
                    int p = ip.getPixel(u, v);
                    if (p > 0) {
                        points_of_single_im.add(Point.create(u, v));
                    }
                }
            }
            points.add(points_of_single_im);
        }
        List<Point> points_X = points.get(0);
        List<Point> points_X_2 = points.get(1);
        RealMatrix transformation = MatrixUtils.createRealMatrix(new double[][]
                    {{0.013, 1.088,18.688},
                    {-1.000,-0.050,127.500}});

        double[][] points_as_array=new double[3][points_X.size()];
        for(int i=0;i<points_X.size();i++)
        {
            points_as_array[0][i]=points_X.get(i).getX();
            points_as_array[1][i]=points_X.get(i).getY();
            points_as_array[2][i]=1;
        }



        AffineTransformator affineTransformator=new AffineTransformator(transformation);
        affineTransformator.transform( MatrixUtils.createRealMatrix(points_as_array));
        ImageProcessor new_im = new ByteProcessor(400, 200);
        RealMatrix results=affineTransformator.getTransformed_points();
        for(int i=0;i<results.getColumnDimension();i++)
        {
                new_im.putPixel((int)results.getEntry(0,i), (int)results.getEntry(1,i), 255);

        }
        (new ImagePlus("abc", new_im)).show();




        triangulator_X = new TriangulationChew(points_X);
        List<Triangle> triangles_X = triangulator_X.getTriangles();
        triangulator_X_2 = new TriangulationChew(points_X_2);

//
//        for (int i = 0; i < triangles_X.size(); i++) {
//            Point[] points_from_X = triangles_X.get(i).getPoints();
//            Point[] points_from_X_prim = triangles_X.get(i).getPoints();
//            double x_x0 = points_from_X[0].getX();
//            double x_y0 = points_from_X[0].getY();
//            double x_x1 = points_from_X[0].getX();
//            double x_y1 = points_from_X[0].getY();
//            double x_x2 = points_from_X[0].getX();
//            double x_y2 = points_from_X[0].getY();
//            double x2_x0 = points_from_X_prim[0].getX();
//            double x2_y0 = points_from_X_prim[0].getY();
//            double x2_x1 = points_from_X_prim[0].getX();
//            double x2_y1 = points_from_X_prim[0].getY();
//            double x2_x2 = points_from_X_prim[0].getX();
//            double x2_y2 = points_from_X_prim[0].getY();
//            RealMatrix A = MatrixUtils.createRealMatrix(new double[][]
//                    {{x_x0,x_y0,1,0,0,0},
//                    {0,0,0,x_x0,x_y0,1},
//                    {x_x1,x_y1,1,0,0,0},
//                    {0,0,0,x_x1,x_y1,1},
//                    {x_x2,x_y2,1,0,0,0},
//                    {0,0,0,x_x2,x_y2,1},
//                    });
//            RealVector b = MatrixUtils.createRealVector(new double[]{ x2_x0, x2_y0, x2_x1, x2_y1, x2_x2, x2_y2);
//            DecompositionSolver s = new SingularValueDecomposition(A).getSolver();
//            RealVector affine_transformation = s.solve(b);
//
//
//
//
//        }
//        }

    }
    }