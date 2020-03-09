package assignment01;

import java.awt.Color;
import java.awt.Point;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;


import ij.IJ;
import ij.ImagePlus;
import ij.plugin.filter.PlugInFilter;
import ij.process.ImageProcessor;

/**
 * This ImageJ plugin collects the foreground points in
 * a binary image, converts to a new color image and
 * draws half of the dots in color.
 *
 * @author WB
 * @version 2020-02-26
 */

//class RandomSampling<T> {
//
//};


public class Assigment_1 implements PlugInFilter {

    public int setup(String arg, ImagePlus im) {
        return PlugInFilter.DOES_8G;
    }
    public static <T> List<T> sample_without_repetition(List<T> arr, int number) {
        List<T> copied_list = new ArrayList<>(arr);
        List<T> final_list = new ArrayList<>();
        Random rand = new Random();
        for (int i = 0; i < number; i++) {
            int index = rand.nextInt(copied_list.size());
            final_list.add(copied_list.get(index));
            copied_list.remove(index);
        }
        return final_list;
    };

    public void run(ImageProcessor ip) {
        final int w = ip.getWidth();
        final int h = ip.getHeight();
        double calculated_x, calculated_y, radius, x1, y1, x2, y2, x3, y3,best_x=0,best_y=0,best_radius=0;
        int ntrials = 2000;
        List<Integer> number_of_points_close_to_circle_in_trials = new ArrayList<>();

        // Collect all image points with pixel values greater than zero:
        List<Point> pntlist = new ArrayList<Point>();
        for (int v = 0; v < h; v++) {
            for (int u = 0; u < w; u++) {
                int p = ip.getPixel(u, v);
                if (p > 0) {
                    pntlist.add(new Point(u, v));
                }
            }
        }
        IJ.log("Found " + pntlist.size() + " foreground points.");

        int max_number_of_points_near_circle=0,number_of_points_close_to_circle=0;
        for (int i = 0; i < ntrials; i++)
        {
            List<Point> selected_points = sample_without_repetition(pntlist, 3);
            x1 = selected_points.get(0).x;
            y1 = selected_points.get(0).y;
            x2 = selected_points.get(1).x;
            y2 = selected_points.get(1).y;
            x3 = selected_points.get(2).x;
            y3 = selected_points.get(2).y;
            double denominator=2 * (x1 * (y2 - y3) - y1 * (x2 - x3) + x2 * y3 - x3 * y2);
            if(denominator!=0)
            {
                calculated_x = ((x1 * x1 + y1 * y1) * (y2 - y3) + (x2 * x2 + y2 * y2) * (y3 - y1) + (x3 * x3 + y3 * y3) * (y1 - y2)) / denominator;
                calculated_y = ((x1 * x1 + y1 * y1) * (x3 - x2) + (x2 * x2 + y2 * y2) * (x1 - x3) + (x3 * x3 + y3 * y3) * (x2 - x1)) / denominator;
                radius = Math.sqrt(Math.pow(calculated_x - x1, 2) + Math.pow(calculated_y - y1, 2));
                number_of_points_close_to_circle = 0;
                for (Point p : pntlist) {
                    if (Math.abs(Math.pow((p.x - calculated_x), 2) + Math.pow((p.y - calculated_y), 2) - Math.pow(radius, 2)) <= 0.5)
                    {
                        number_of_points_close_to_circle++;
                    }
                }
                if(number_of_points_close_to_circle>=max_number_of_points_near_circle)
                {
                    max_number_of_points_near_circle=number_of_points_close_to_circle;
                    best_x=calculated_x;
                    best_y=calculated_y;
                    best_radius=radius;
                }
                number_of_points_close_to_circle_in_trials.add(number_of_points_close_to_circle);
            }
        }

        ImageProcessor cp = ip.convertToColorProcessor();
        cp.setColor(Color.red);
        cp.drawOval((int)(best_x-best_radius), (int)(best_y-best_radius), (int)(2*best_radius), (int)(2*best_radius));

//        for (int i = 0; i < 100; i++) {
//            double phi = rg.nextDouble() * 2 * Math.PI;	// uniformly distributed angle 0,...,2pi
//            int u = (int) Math.round(best_x + best_radius * Math.cos(phi) );
//            int v = (int) Math.round(best_y + best_radius * Math.sin(phi) );
//            cp.putPixel(u, v, 255);
//        }
////
//        // Copy 'ip' to a new color image and redraw some of the dots in red:
//        ImageProcessor cp = ip.convertToColorProcessor();
//        cp.setColor(Color.red);
//        for (Point p : pntlist) {
//            //if (p.y <= h / 2)
//                cp.drawDot(p.x, p.y);
//        }
//
//        // Just for show, draw a blue circle somewhere:
//        cp.setColor(Color.blue);
//        cp.setLineWidth(1);
//        cp.drawOval(35, 200, 75, 75);
//
//        // Display the newly created image:
        showImage(cp, "colored dots");
    }

    void showImage(ImageProcessor ip, String title) {
        (new ImagePlus(title, ip)).show();
    }

}
