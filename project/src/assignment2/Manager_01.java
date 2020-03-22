package assignment2;

import ij.IJ;
import ij.ImagePlus;
import ij.io.Opener;
import ij.plugin.filter.PlugInFilter;
import ij.process.ByteProcessor;
import ij.process.ImageProcessor;
import ij.ImageStack;
import imagingbook.pub.geometry.basic.Point;
import imagingbook.pub.geometry.delaunay.chew.TriangulationChew;
import imagingbook.pub.geometry.delaunay.Triangle;

import java.util.ArrayList;

import java.util.List;


import org.apache.commons.math3.linear.*;
import org.paukov.combinatorics3.Generator;

class AffineTransformator_ {
    private RealMatrix transformations;
    private RealMatrix transformed_points;

    public AffineTransformator_(RealMatrix transformations) {
        this.transformations = transformations;
    }


    public void setTransformations(RealMatrix transformations) {
        this.transformations = transformations;
    }

    public RealMatrix getTransformed_points() {
        return transformed_points;
    }

    public void transform(RealMatrix points) {
        this.transformed_points = this.transformations.multiply(points);
    }

    public double calculate_error(RealMatrix points) {

        double error = 0;
        double min_distance;
        double distance;
        for (int i = 0; i < this.transformed_points.getColumnDimension(); i++) {
            min_distance = 1000000;
            for (int j = 0; j < points.getColumnDimension(); j++)
            {
                distance = Math.sqrt(Math.pow(this.transformed_points.getEntry(0, i) - points.getEntry(0, j), 2) + Math.pow(this.transformed_points.getEntry(1, i) - points.getEntry(1, j), 2));
                if (distance < min_distance)
                {
                    min_distance = distance;
                }
            }
            error += min_distance;
        }
        return error;
    }


}

public class Manager_01 implements PlugInFilter {
    private TriangulationChew triangulator_X;
    private TriangulationChew triangulator_X_2;
    ImagePlus im = null;


    public int setup(String arg, ImagePlus im) {

        this.im = im;
        return DOES_8G + STACK_REQUIRED;
    }
    List<Point> getPointsFromProcessor(ImageProcessor improc)
    {
        List<Point> points_of_single_im = new ArrayList<Point>();
        final int w = improc.getWidth();
        final int h = improc.getHeight();
        for (int v = 0; v < h; v++) {
            for (int u = 0; u < w; u++) {
                int p = improc.getPixel(u, v);
                if (p > 0) {
                    points_of_single_im.add(Point.create(u, v));

                }
            }
        }
        return points_of_single_im;
    }
    public void run(ImageProcessor ip) {

        ImageStack stack = im.getStack();

        // Collect all image points with pixel values greater than zero:
//        List<List<Point>> points = new ArrayList<List<Point>>();
//        for (int k = 1; k <= 2; k++) {    // NOTE: k starts from 1!!
//            ImageProcessor ipk = stack.getProcessor(k);
//            List<Point> points_of_single_im = new ArrayList<Point>();
//            IJ.log("Reading stack nr"+k);
//
//            final int w = ipk.getWidth();
//            final int h = ipk.getHeight();
//            for (int v = 0; v < h; v++) {
//                for (int u = 0; u < w; u++) {
//                    int p = ip.getPixel(u, v);
//                    if (p > 0) {
//                        points_of_single_im.add(Point.create(u, v));
//                        IJ.log("Point of stack"+k+"with coordinate"+u+" "+v);
//
//                    }
//                }
//            }
//            points.add(points_of_single_im);
//        }



        Opener opener = new Opener();
        String imageFilePath_X = "C:/Users/blaze/Desktop/FHOO/assignment_02_computer_vision/test-point-set1.png";
        String imageFilePath_X_prim = "C:/Users/blaze/Desktop/FHOO/assignment_02_computer_vision/test-point-set2.png";
        ImagePlus imp_X = opener.openImage(imageFilePath_X);
        ImagePlus imp_X_prim = opener.openImage(imageFilePath_X_prim);
        ImageProcessor imp_X_processor = imp_X.getProcessor();
        ImageProcessor imp_X_prim_processor = imp_X_prim.getProcessor();

        List<Point> points_X = getPointsFromProcessor(imp_X_processor);
        List<Point> points_X_prim = getPointsFromProcessor(imp_X_prim_processor);

//        List<Point> points_X = points.get(0);
//        List<Point> points_X_prim = points.get(1);

        RealMatrix transformation = MatrixUtils.createRealMatrix(new double[][]
                {{0.013, 1.088, 18.688},
                        {-1.000, -0.050, 127.500}});

        double[][] points_as_array_X = new double[3][points_X.size()];
        for (int i = 0; i < points_X.size(); i++) {
            points_as_array_X[0][i] = points_X.get(i).getX();
            points_as_array_X[1][i] = points_X.get(i).getY();
            points_as_array_X[2][i] = 1;
        }

        double[][] points_as_array_X_prim = new double[3][points_X_prim.size()];
        for (int i = 0; i < points_X_prim.size(); i++) {
            points_as_array_X_prim[0][i] = points_X_prim.get(i).getX();
            points_as_array_X_prim[1][i] = points_X_prim.get(i).getY();
            points_as_array_X_prim[2][i] = 1;
        }


        AffineTransformator_ affineTransformator = new AffineTransformator_(transformation);
        affineTransformator.transform(MatrixUtils.createRealMatrix(points_as_array_X));
        for (int i = 0; i < points_X_prim.size(); i++) {
//            IJ.log("X " + String.valueOf(points_X_prim.get(i).getX()));
//            IJ.log("Y " + String.valueOf(points_X_prim.get(i).getY()));
            IJ.log("X' " + String.valueOf(points_as_array_X_prim[0][i]));
            IJ.log("Y' " + String.valueOf(points_as_array_X_prim[1][i]));

        }

        RealMatrix results = affineTransformator.getTransformed_points();

        for (int i = 0; i < results.getColumnDimension(); i++) {
            IJ.log("transformed X" + String.valueOf((int)results.getEntry(0, i)));
            IJ.log("transformed Y" + String.valueOf((int)results.getEntry(1, i)));

        }
        IJ.log("default affine ");
        for (int i = 0; i < transformation.getRowDimension(); i++) {
            for(int j=0;j<transformation.getColumnDimension();j++)
                IJ.log(String.valueOf(transformation.getEntry(i,j)));
        }

        IJ.log("error " + String.valueOf(affineTransformator.calculate_error(MatrixUtils.createRealMatrix(points_as_array_X_prim).getSubMatrix(0, 1, 0, points_X_prim.size() - 1))));


        ImageProcessor new_im = new ByteProcessor(400, 200);

        for (int i = 0; i < results.getColumnDimension(); i++) {
            new_im.putPixel((int) results.getEntry(0, i), (int) results.getEntry(1, i), 255);

        }
        (new ImagePlus("transformed points", new_im)).show();


        triangulator_X = new TriangulationChew(points_X);
        List<Triangle> triangles_X = triangulator_X.getTriangles();
        triangulator_X_2 = new TriangulationChew(points_X_prim);
        List<Triangle> triangles_X_prim = triangulator_X_2.getTriangles();
        double best_triangulation_error = 1000000, current_triangulation_error = 0;
        RealVector bestAffine = new ArrayRealVector();
        RealMatrix bestProjection = new Array2DRowRealMatrix();

        for (int iter_triangles_X = 0; iter_triangles_X < triangles_X.size(); iter_triangles_X++) {
            Point[] points_from_X = triangles_X.get(iter_triangles_X).getPoints();
            double x_x0 = points_from_X[0].getX();
            double x_y0 = points_from_X[0].getY();
            double x_x1 = points_from_X[1].getX();
            double x_y1 = points_from_X[1].getY();
            double x_x2 = points_from_X[2].getX();
            double x_y2 = points_from_X[2].getY();
            for (int iter_triangles_y = 0; iter_triangles_y < triangles_X_prim.size(); iter_triangles_y++) {

                Point[] points_from_X_prim = triangles_X_prim.get(iter_triangles_y).getPoints();
                for(List<Point> permutator: Generator.permutation(points_from_X_prim).simple()) {

                    double x2_x0 = permutator.get(0).getX();
                    double x2_y0 = permutator.get(0).getY();
                    double x2_x1 = permutator.get(1).getX();
                    double x2_y1 = permutator.get(1).getY();
                    double x2_x2 = permutator.get(2).getX();
                    double x2_y2 = permutator.get(2).getY();

                    RealMatrix A = MatrixUtils.createRealMatrix(new double[][]
                            {{x_x0, x_y0, 1, 0, 0, 0},
                                    {0, 0, 0, x_x0, x_y0, 1},
                                    {x_x1, x_y1, 1, 0, 0, 0},
                                    {0, 0, 0, x_x1, x_y1, 1},
                                    {x_x2, x_y2, 1, 0, 0, 0},
                                    {0, 0, 0, x_x2, x_y2, 1},
                            });

                    RealVector b = MatrixUtils.createRealVector(new double[]{x2_x0, x2_y0, x2_x1, x2_y1, x2_x2, x2_y2});
                    DecompositionSolver s = new SingularValueDecomposition(A).getSolver();
                    RealVector affine_transformation = s.solve(b);
                    affineTransformator.setTransformations(MatrixUtils.createRealMatrix(new double[][]{{affine_transformation.getEntry(0),
                            affine_transformation.getEntry(1),
                            affine_transformation.getEntry(2)
                    },
                            {affine_transformation.getEntry(3),
                                    affine_transformation.getEntry(4),
                                    affine_transformation.getEntry(5)
                            }
                    }));
                    affineTransformator.transform(MatrixUtils.createRealMatrix(points_as_array_X));
                    current_triangulation_error = affineTransformator.calculate_error(MatrixUtils.createRealMatrix(points_as_array_X_prim).getSubMatrix(0, 1, 0, points_X_prim.size() - 1));
                    if (current_triangulation_error < best_triangulation_error) {
                        best_triangulation_error = current_triangulation_error;
                        bestAffine = affine_transformation;
                        bestProjection=affineTransformator.getTransformed_points();
                    }
                }
            }
        }

        IJ.log("best affine from triangulation:");

        for (int i = 0; i < bestAffine.getDimension(); i++) {
            IJ.log(String.valueOf(bestAffine.getEntry(i)));
        }




    }
}