package fem.solver.gui;

import fem.solver.FEM;
import fem.solver.Solution;
import javafx.application.Application;
import javafx.scene.Scene;
import javafx.scene.chart.LineChart;
import javafx.scene.chart.NumberAxis;
import javafx.scene.chart.XYChart;
import javafx.scene.layout.BorderPane;
import javafx.stage.Stage;

public class App extends Application {

    private FEM femSolver;
    private LineChart<Number,Number> plot;
    private BorderPane borderPane;
    private final int n = 10;
    private final int domainStart = 0;
    private final int domainEnd = 3;


    public static void main(String[] args) {
        launch(args);
    }

    @Override
    public void init() throws Exception {
        super.init();
        this.femSolver = new FEM(n,domainStart,domainEnd);

        final NumberAxis xAxis = new NumberAxis();
        final NumberAxis yAxis = new NumberAxis();
        this.plot = new LineChart<Number,Number>(xAxis,yAxis);

        borderPane = new BorderPane();
        this.calculate();
        borderPane.setCenter(plot);
    }

    @Override
    public void start(Stage stage) throws Exception {


        Scene scene = new Scene(borderPane,1280,720);

        stage.setScene(scene);
        stage.show();

    }

    private void calculate(){
        Solution solution = this.femSolver.solve();
        XYChart.Series<Number,Number> series = new XYChart.Series<>();

        for (int i=0;i< solution.result().length;i++){
            double x = (solution.domainEnd() - solution.domainStart()) * i / (solution.result().length-1);
            series.getData().add(new XYChart.Data<>(x,solution.result()[i]));
        }
        if(this.plot.getData().size() > 1){
            this.plot.getData().remove(1);
        }
        this.plot.getData().add(series);

    }
}
