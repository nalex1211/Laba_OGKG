package org.example.laba_ogkg;

import javafx.application.Application;
import javafx.scene.Scene;
import javafx.scene.control.Button;
import javafx.scene.control.Label;
import javafx.scene.control.TextField;
import javafx.scene.input.MouseButton;
import javafx.scene.layout.BorderPane;
import javafx.scene.layout.HBox;
import javafx.scene.layout.Pane;
import javafx.scene.paint.Color;
import javafx.scene.shape.Circle;
import javafx.scene.shape.Ellipse;
import javafx.scene.shape.Line;
import javafx.stage.Stage;

import org.apache.commons.math3.linear.*;
import org.apache.commons.math3.random.RandomDataGenerator;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;

public class EllipseApplication extends Application {
    private List<double[]> points = new ArrayList<>();
    private Pane drawingPane = new Pane();
    private double lastX, lastY;
    private double scale = 1.0;
    private double width = 1920;
    private double height = 1080;
    private List<javafx.scene.shape.Shape> shapes = new ArrayList<>();

    @Override
    public void start(Stage primaryStage) {
        BorderPane root = new BorderPane();

        HBox inputBar = new HBox(10);
        Label label = new Label("Кількість Точок:");
        TextField textField = new TextField();
        Button btnGeneratePoints = new Button("Згенерувати Точки");
        btnGeneratePoints.setOnAction(e -> generatePoints(Integer.parseInt(textField.getText()), drawingPane.getWidth(), drawingPane.getHeight()));

        inputBar.getChildren().addAll(label, textField, btnGeneratePoints);
        root.setTop(inputBar);

        HBox buttonBar = new HBox(10);
        Button btnGenerateEllipse = new Button("Згенерувати Еліпс");
        btnGenerateEllipse.setOnAction(e -> {
            clearShapes();
            List<double[]> outerShell = getConvexHull(points.toArray(new double[0][]));
            if (!outerShell.isEmpty()) {
                drawConvexHull(outerShell);
                double[][] shellArray = outerShell.toArray(new double[0][]);
                RealMatrix results = mvee(shellArray, 0.001);
                RealMatrix A = results.getSubMatrix(0, 1, 0, 1);
                double[] centroid = results.getColumn(2);
                drawEllipse(A, centroid);
            }
        });

        Button btnClear = new Button("Очистити Полотно");
        btnClear.setOnAction(e -> clearCanvas());

        buttonBar.getChildren().addAll(btnGenerateEllipse, btnClear);
        root.setBottom(buttonBar);

        root.setCenter(drawingPane);
        enableCanvasDraggingAndScaling(root);

        drawingPane.setOnMouseClicked(e -> {
            if (e.getButton() == MouseButton.PRIMARY) {
                double x = e.getX();
                double y = e.getY();
                points.add(new double[]{x, y});
                Circle circle = new Circle(x, y, 3, Color.RED);
                drawingPane.getChildren().add(circle);
            }
        });

        Scene scene = new Scene(root, width, height);
        primaryStage.setTitle("Мінімальний Охоплюючий Еліпс");
        primaryStage.setScene(scene);
        primaryStage.show();
    }

    private void clearCanvas() {
        drawingPane.getChildren().clear();
        points.clear();
        shapes.clear();
    }

    private void clearShapes() {
        drawingPane.getChildren().removeAll(shapes);
        shapes.clear();
    }

    private void generatePoints(int numPoints, double width, double height) {
        clearCanvas();
        RandomDataGenerator randomDataGenerator = new RandomDataGenerator();
        double meanX = width / 2;
        double meanY = height / 2;
        double stdDevX = width / 6;
        double stdDevY = height / 6;

        for (int i = 0; i < numPoints; i++) {
            double x = randomDataGenerator.nextGaussian(meanX, stdDevX);
            double y = randomDataGenerator.nextGaussian(meanY, stdDevY);

            x = Math.max(0, Math.min(width, x));
            y = Math.max(0, Math.min(height, y));

            points.add(new double[]{x, y});
            Circle circle = new Circle(x, y, 3, Color.RED);
            drawingPane.getChildren().add(circle);
        }
    }

    private void drawEllipse(RealMatrix A, double[] centroid) {
        SingularValueDecomposition svd = new SingularValueDecomposition(A);
        double[] singularValues = svd.getSingularValues();
        RealMatrix uMatrix = svd.getU();

        double rx = 1.0 / Math.sqrt(singularValues[0]);
        double ry = 1.0 / Math.sqrt(singularValues[1]);

        double angle = Math.atan2(uMatrix.getEntry(1, 0), uMatrix.getEntry(0, 0));

        Ellipse ellipse = new Ellipse(centroid[0], centroid[1], rx, ry);
        ellipse.setRotate(Math.toDegrees(angle));
        ellipse.setStroke(Color.BLUE);
        ellipse.setFill(Color.TRANSPARENT);
        drawingPane.getChildren().add(ellipse);
        shapes.add(ellipse);
    }

    private void enableCanvasDraggingAndScaling(BorderPane root) {
        root.setOnMousePressed(e -> {
            if (e.getButton() == MouseButton.SECONDARY) {
                lastX = e.getSceneX();
                lastY = e.getSceneY();
            }
        });

        root.setOnMouseDragged(e -> {
            if (e.getButton() == MouseButton.SECONDARY) {
                double deltaX = e.getSceneX() - lastX;
                double deltaY = e.getSceneY() - lastY;
                drawingPane.setTranslateX(drawingPane.getTranslateX() + deltaX);
                drawingPane.setTranslateY(drawingPane.getTranslateY() + deltaY);
                lastX = e.getSceneX();
                lastY = e.getSceneY();
            }
        });

        root.setOnScroll(e -> {
            double deltaY = e.getDeltaY();
            scale += deltaY * 0.001;

            drawingPane.setScaleX(scale);
            drawingPane.setScaleY(scale);

        });
    }

    public static RealMatrix mvee(double[][] points, double tol) {
        List<double[]> hull = getConvexHull(points);
        int N = hull.size();
        int d = points[0].length;

        RealMatrix Q = MatrixUtils.createRealMatrix(d + 1, N);
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < d; j++) {
                Q.setEntry(j, i, hull.get(i)[j]);
            }
            Q.setEntry(d, i, 1);
        }

        RealMatrix X = null;
        RealVector u = new ArrayRealVector(N, 1.0 / N);
        double err = tol + 1.0;
        int maxIterations = 1000;
        int iteration = 0;

        while (err > tol && iteration < maxIterations) {
            X = Q.multiply(new DiagonalMatrix(u.toArray())).multiply(Q.transpose());
            RealMatrix M_inv = new LUDecomposition(X).getSolver().getInverse();
            RealMatrix M = Q.transpose().multiply(M_inv).multiply(Q);
            double[] diagM = getDiagonal(M);
            int jdx = getMaxIndex(diagM);
            double max = diagM[jdx];
            double step_size = (max - d - 1.0) / ((d + 1) * (max - 1.0));
            RealVector new_u = u.mapMultiply(1 - step_size);
            new_u.setEntry(jdx, new_u.getEntry(jdx) + step_size);
            err = new_u.subtract(u).getNorm();
            u = new_u;
            iteration++;
        }

        RealVector c = Q.operate(u);
        RealMatrix outerProduct = c.outerProduct(c);
        RealMatrix subtracted = X.subtract(outerProduct);
        RealMatrix reducedSubtracted = subtracted.getSubMatrix(0, d - 1, 0, d - 1);
        RealMatrix A = new LUDecomposition(reducedSubtracted).getSolver().getInverse().scalarMultiply(1.0 / d);

        double[] centroid = Arrays.copyOf(c.toArray(), d);
        return MatrixUtils.createRealMatrix(new double[][]{
                {A.getEntry(0, 0), A.getEntry(0, 1), centroid[0]},
                {A.getEntry(1, 0), A.getEntry(1, 1), centroid[1]}
        });
    }

    private static List<double[]> getConvexHull(double[][] points) {
        List<double[]> sortedPoints = Arrays.asList(points.clone());
        sortedPoints.sort(Comparator.comparingDouble(p -> p[0]));

        List<double[]> lowerHull = new ArrayList<>();
        for (double[] p : sortedPoints) {
            while (lowerHull.size() >= 2 && cross(lowerHull.get(lowerHull.size() - 2), lowerHull.get(lowerHull.size() - 1), p) <= 0) {
                lowerHull.remove(lowerHull.size() - 1);
            }
            lowerHull.add(p);
        }

        List<double[]> upperHull = new ArrayList<>();
        for (int i = sortedPoints.size() - 1; i >= 0; i--) {
            double[] p = sortedPoints.get(i);
            while (upperHull.size() >= 2 && cross(upperHull.get(upperHull.size() - 2), upperHull.get(upperHull.size() - 1), p) <= 0) {
                upperHull.remove(upperHull.size() - 1);
            }
            upperHull.add(p);
        }

        lowerHull.remove(lowerHull.size() - 1);
        upperHull.remove(upperHull.size() - 1);
        lowerHull.addAll(upperHull);

        return lowerHull;
    }

    private static double[] getDiagonal(RealMatrix matrix) {
        int minDim = Math.min(matrix.getRowDimension(), matrix.getColumnDimension());
        double[] diag = new double[minDim];
        for (int i = 0; i < minDim; i++) {
            diag[i] = matrix.getEntry(i, i);
        }
        return diag;
    }

    private static int getMaxIndex(double[] array) {
        int maxIndex = 0;
        double maxValue = array[0];
        for (int i = 1; i < array.length; i++) {
            if (array[i] > maxValue) {
                maxValue = array[i];
                maxIndex = i;
            }
        }
        return maxIndex;
    }

    private static double cross(double[] o, double[] a, double[] b) {
        return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0]);
    }

    private void drawConvexHull(List<double[]> hull) {
        if (hull.size() < 2) return;

        for (int i = 0; i < hull.size(); i++) {
            double[] p1 = hull.get(i);
            double[] p2 = hull.get((i + 1) % hull.size());
            Line line = new Line(p1[0], p1[1], p2[0], p2[1]);
            line.setStroke(Color.GREEN);
            drawingPane.getChildren().add(line);
            shapes.add(line);
        }
    }

    public static void main(String[] args) {
        launch(args);
    }
}
