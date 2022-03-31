package fem.solver;

import java.util.Arrays;

public record Solution(double domainStart, double domainEnd, double[] result) {

    @Override
    public String toString() {
        return Arrays.toString(result);
    }
}

