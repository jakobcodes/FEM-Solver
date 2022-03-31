package fem.solver;

import org.apache.commons.math3.analysis.integration.IterativeLegendreGaussIntegrator;
import org.apache.commons.math3.analysis.integration.UnivariateIntegrator;
import org.apache.commons.math3.linear.*;

import java.util.Arrays;

public class FEM {

    public static final double G = 6.6743e-11;
    private final int n;
    private final int domainStart;
    private final int domainEnd;
    private final double h;
    private final UnivariateIntegrator integrator;


    public FEM(int n,int domainStart, int domainEnd) {
        this.n = n;
        this.domainStart = domainStart;
        this.domainEnd = domainEnd;
        this.h = (domainEnd - domainStart)/(double) n;

        // INTEGRATOR
        int integrationPoints = 5;
        this.integrator = new IterativeLegendreGaussIntegrator(
                integrationPoints,
                1e-6,
                1e-6
        );
    }

    public Solution solve(){
        RealMatrix bMatrix = new Array2DRowRealMatrix(n+1,n+1);

        for (int i=1; i < n;i++){
            for (int j=i-1;j <= i;j++){
                if(j == 0 ) continue;
                if(Math.abs(i - j) > 1) continue;

                double x0,x2;
                if(i == j) {
                    x0 = (i-1) *this.h;
                    x2 = (i+1) *this.h;
                }else{
                    x0 = j *this.h;
                    x2 = i *this.h;
                }
                int finalJ = j;
                int finalI = i;
                double integral = this.integrator.integrate(
                        Integer.MAX_VALUE,
                        x -> de_i(finalI,x) * de_i(finalJ,x),
                        x0,
                        x2
                );
                bMatrix.setEntry(i,j,-integral);
                bMatrix.setEntry(j,i,-integral);
            }
        }
        bMatrix.setEntry(0,0,1);
        bMatrix.setEntry(n,n,1);

        // build L vector
        RealVector lVector = new ArrayRealVector(n+1,0);

        for(int i=1;i < n;i++){
            double x0 = (i-1)*this.h;
            double x2 = (i+1)*this.h;
            int finalI = i;
            double integral1 = this.integrator.integrate(
                    Integer.MAX_VALUE,
                    x -> 4*Math.PI * G * p(x) * e_i(finalI,x),
                    x0,
                    x2 - this.h
            );
            double integral2 = this.integrator.integrate(
                    Integer.MAX_VALUE,
                    x -> 4*Math.PI * G * p(x) * e_i(finalI,x),
                    x0+this.h,
                    x2
            );
            double integral3 = this.integrator.integrate(
                    Integer.MAX_VALUE,
                    x -> deshiftU(x)* de_i(finalI, x),
                    x0,
                    x2 - this.h
            );
            double integral4 = this.integrator.integrate(
                    Integer.MAX_VALUE,
                    x -> deshiftU(x)* de_i(finalI, x),
                    x0 + this.h,
                    x2
            );
            lVector.setEntry(i,integral1 + integral2 + integral3 + integral4);
        }
        // calculate wVector
        RealVector wVector = new LUDecomposition(bMatrix).getSolver().solve(lVector);

        // calculate result Vector
        double[] result = Arrays.copyOf(wVector.toArray(), n+1);
        double[] xArray = new double[n+1];
        for(int i=0;i<n+1;i++){
            xArray[i] = i * this.h;
        }

        result[n] = 0;
        for(int i=0;i<n+1;i++){
            result[i] += shiftU(xArray[i]);
            System.out.println(result[i]);
        }

        System.out.println(bMatrix);
        System.out.println(lVector);
        return new Solution(domainStart,domainEnd,result);

    }

    public double p(double x){
        if (x > 1 && x <= 2) return 1;
        else return 0;
    }

    public double e_i(int i, double x){
        double x1 = h * (i-1);
        double x2 = h * (i);
        double x3 = h * (i+1);

        if(x > x1 && x <= x2) return (x-x1)/this.h;
        else if (x > x2 && x < x3) return (x3 - x)/this.h;
        else return 0;
    }
    public double de_i(int i, double x){
        double x1 = h * (i-1);
        double x2 = h * (i);
        double x3 = h * (i+1);

        if(x > x1 && x <= x2) return 1/this.h;
        else if (x > x2 && x < x3) return -1/this.h;
        else return 0;
    }
    public double shiftU(double x){
        return (-x/3) + 5;
    }
    public double deshiftU(double x){
        return (double) (-1/3) + 5;
    }
}
