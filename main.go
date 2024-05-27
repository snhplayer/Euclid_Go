package main

import (
    "fmt"
    "math/big"
    "math/rand"
    "strings"
    "time"
    "gonum.org/v1/plot"
    "gonum.org/v1/plot/plotter"
    "gonum.org/v1/plot/vg"
)

// polyRing represents a polynomial ring over rational numbers
type polyRing struct {
    coeff []*big.Rat
}

// newPolyRing creates a new polynomial from the given coefficients
func newPolyRing(coeffs []*big.Rat) *polyRing {
    return &polyRing{coeff: coeffs}
}

// deg returns the degree of the polynomial
func (p *polyRing) deg() int {
    for i := len(p.coeff) - 1; i >= 0; i-- {
        if p.coeff[i].Sign() != 0 {
            return i
        }
    }
    return 0
}

// isZero checks if the polynomial is zero
func (p *polyRing) isZero() bool {
    for _, c := range p.coeff {
        if c.Sign() != 0 {
            return false
        }
    }
    return true
}

func (p *polyRing) String() string {
    var b strings.Builder
    for i := len(p.coeff) - 1; i >= 0; i-- {
        if p.coeff[i].Sign() != 0 {
            if b.Len() > 0 && p.coeff[i].Sign() > 0 {
                b.WriteString(" + ")
            } else if p.coeff[i].Sign() < 0 {
                b.WriteString(" - ")
            }
            if absRat(p.coeff[i]).Cmp(big.NewRat(1, 1)) != 0 || i == 0 {
                b.WriteString(absRat(p.coeff[i]).String())
                if i > 0 {
                    b.WriteString("*")
                }
            }
            if i > 0 {
                b.WriteString("x")
                if i > 1 {
                    b.WriteString("^" + fmt.Sprint(i))
                }
            }
        }
    }
    return b.String()
}

func absRat(r *big.Rat) *big.Rat {
    if r.Sign() < 0 {
        return new(big.Rat).Neg(r)
    }
    return r
}

// add adds two polynomials
func (p *polyRing) add(q *polyRing) *polyRing {
    maxDeg := max(p.deg(), q.deg())
    result := make([]*big.Rat, maxDeg+1)
    for i := 0; i <= maxDeg; i++ {
        result[i] = new(big.Rat)
        if i <= p.deg() {
            result[i].Add(result[i], p.coeff[i])
        }
        if i <= q.deg() {
            result[i].Add(result[i], q.coeff[i])
        }
    }
    return newPolyRing(result)
}

// sub subtracts two polynomials
func (p *polyRing) sub(q *polyRing) *polyRing {
    maxDeg := max(p.deg(), q.deg())
    result := make([]*big.Rat, maxDeg+1)
    for i := 0; i <= maxDeg; i++ {
        result[i] = new(big.Rat)
        if i <= p.deg() {
            result[i].Add(result[i], p.coeff[i])
        }
        if i <= q.deg() {
            result[i].Sub(result[i], q.coeff[i])
        }
    }
    return newPolyRing(result)
}

// mul multiplies two polynomials
func (p *polyRing) mul(q *polyRing) *polyRing {
    result := make([]*big.Rat, p.deg()+q.deg()+1)
    for i := range result {
        result[i] = new(big.Rat)
    }
    for i := range p.coeff {
        for j := range q.coeff {
            temp := new(big.Rat).Mul(p.coeff[i], q.coeff[j])
            result[i+j].Add(result[i+j], temp)
        }
    }
    return newPolyRing(result)
}

func (p *polyRing) div(q *polyRing) (*polyRing, *polyRing) {
    if q.isZero() {
        panic("division by zero")
    }

    pDeg, qDeg := p.deg(), q.deg()
    if pDeg < qDeg {
        // If the degree of p is less than the degree of q, return quotient as 0 and p as the remainder
        return newPolyRing([]*big.Rat{new(big.Rat)}), newPolyRing(p.coeff)
    }

    quotient := make([]*big.Rat, pDeg-qDeg+1)
    remainder := make([]*big.Rat, pDeg+1) // Ensure this matches the degree of p

    for i := range quotient {
        quotient[i] = new(big.Rat)
    }
    for i := range remainder {
        remainder[i] = new(big.Rat).Set(p.coeff[i])
    }

    for pDeg >= qDeg {
        leadCoeff := new(big.Rat).Quo(remainder[pDeg], q.coeff[qDeg])
        quotient[pDeg-qDeg] = leadCoeff

        for i := range q.coeff {
            temp := new(big.Rat).Mul(leadCoeff, q.coeff[i])
            if pDeg-qDeg+i < len(remainder) {
                remainder[pDeg-qDeg+i].Sub(remainder[pDeg-qDeg+i], temp)
            }
        }

        for pDeg >= 0 && remainder[pDeg].Sign() == 0 {
            pDeg--
        }
    }

    // Ensure the remainder slice is correctly sliced to match the actual degree
    return newPolyRing(quotient), newPolyRing(remainder[:pDeg+1])
}





// extendedEuclideanPoly implements the extended Euclidean algorithm for polynomials
func extendedEuclideanPoly(f, g *polyRing) (*polyRing, *polyRing, *polyRing) {
    s0 := newPolyRing([]*big.Rat{big.NewRat(1, 1)})
    s1 := newPolyRing([]*big.Rat{new(big.Rat)})
    t0 := newPolyRing([]*big.Rat{new(big.Rat)})
    t1 := newPolyRing([]*big.Rat{big.NewRat(1, 1)})

    for !g.isZero() {
        q, r := f.div(g)
        f, g = g, r
        s0, s1 = s1, s0.sub(q.mul(s1))
        t0, t1 = t1, t0.sub(q.mul(t1))
    }

    return f, s0, t0
}

func max(a, b int) int {
    if a > b {
        return a
    }
    return b
}

func generateRandomPolynomial(degree int) *polyRing {
    coeffs := make([]*big.Rat, degree+1)
    for i := 0; i <= degree; i++ {
        coeff := rand.Intn(11) - 5 // Random coefficient between -5 and 5
        coeffs[i] = big.NewRat(int64(coeff), 1)
    }
    return newPolyRing(coeffs)
}

func colorize(text, color string) string {
    return fmt.Sprintf("%s%s%s", color, text, "\033[0m")
}

func testExtendedEuclidean(numTests int) {
    for i := 0; i < numTests; i++ {
        degreeF := rand.Intn(5) + 1 // Random degree between 1 and 5
        degreeG := rand.Intn(5) + 1 // Random degree between 1 and 5

        f := generateRandomPolynomial(degreeF)
        g := generateRandomPolynomial(degreeG)

        // Ensure g is not zero
        for g.isZero() {
            g = generateRandomPolynomial(degreeG)
        }

        startTime := time.Now()

        // Perform extended Euclidean algorithm
        gcd, s, t := extendedEuclideanPoly(f, g)

        endTime := time.Now()
        totalTime := endTime.Sub(startTime)

        // Print results
        fmt.Printf("\n%s %d\n", colorize("Test", "\033[1;34m"), i+1)
        fmt.Printf("%s %v\n", colorize("f(x):", "\033[1;32m"), f)
        fmt.Printf("%s %v\n", colorize("g(x):", "\033[1;32m"), g)
        fmt.Printf("%s %v\n", colorize("GCD:", "\033[1;33m"), gcd)
        fmt.Printf("%s %v\n", colorize("s(x):", "\033[1;36m"), s)
        fmt.Printf("%s %v\n", colorize("t(x):", "\033[1;36m"), t)
        fmt.Printf("%s %.6f seconds\n", colorize("Execution time:", "\033[1;35m"), totalTime.Seconds())
    }
}

func testExtendedEuclideanLength(maxLength int) {
    points := make(plotter.XYs, maxLength)
    var totalTime time.Duration

    for i := 1; i <= maxLength; i++ {
        f := generateRandomPolynomial(i)
        g := generateRandomPolynomial(i)

        startTime := time.Now()
        extendedEuclideanPoly(f, g)
        endTime := time.Now()
        totalTime += endTime.Sub(startTime)

        points[i-1].X = float64(i)
        points[i-1].Y = totalTime.Seconds()
    }

    fmt.Printf("%s %.6f seconds\n", colorize("Total execution time:", "\033[1;35m"), totalTime.Seconds())

    p := plot.New()
    p.Title.Text = "Polynomial Length vs. Execution Time"
    p.X.Label.Text = "Polynomial Length"
    p.Y.Label.Text = "Execution Time (seconds)"

    line, err := plotter.NewLine(points)
    if err != nil {
        panic(err)
    }
    p.Add(line)

    // Save the plot to a PNG file.
    if err := p.Save(6*vg.Inch, 4*vg.Inch, "plot.png"); err != nil {
        panic(err)
    }
}



func main() {
    rand.Seed(time.Now().UnixNano())

    // Input coefficients of the first polynomial
    fmt.Print("Enter the degree of the first polynomial: ")
    var degreeF int
    fmt.Scanln(&degreeF)

    coeffsF := make([]*big.Rat, degreeF+1)
    for i := degreeF; i >= 0; i-- {
        fmt.Printf("Enter the coefficient for x^%d: ", i)
        var coeff int64
        fmt.Scanln(&coeff)
        coeffsF[i] = big.NewRat(coeff, 1)
    }
    f := newPolyRing(coeffsF)

    // Input coefficients of the second polynomial
    fmt.Print("Enter the degree of the second polynomial: ")
    var degreeG int
    fmt.Scanln(&degreeG)

    coeffsG := make([]*big.Rat, degreeG+1)
    for i := degreeG; i >= 0; i-- {
        fmt.Printf("Enter the coefficient for x^%d: ", i)
        var coeff int64
        fmt.Scanln(&coeff)
        coeffsG[i] = big.NewRat(coeff, 1)
    }
    g := newPolyRing(coeffsG)

    // Start timing
    startTime := time.Now()

    // Perform extended Euclidean algorithm
    gcd, s, t := extendedEuclideanPoly(f, g)

    // End timing
    endTime := time.Now()
    totalTime := endTime.Sub(startTime)

    // Print results
    fmt.Printf("\n%s %v\n", colorize("GCD of the two polynomials:", "\033[1;33m"), gcd)
    fmt.Printf("%s %v\n", colorize("U(x):", "\033[1;36m"), s)
    fmt.Printf("%s %v\n", colorize("V(x):", "\033[1;36m"), t)
    fmt.Printf("%s %.6f seconds\n", colorize("Execution time:", "\033[1;35m"), totalTime.Seconds())

    // Run tests
    fmt.Print("\nEnter the number of random tests to run: ")
    var numTests int
    fmt.Scanln(&numTests)
    testExtendedEuclidean(numTests)

    fmt.Print("\nEnter the length of random polynoms to test: ")
    var numTestsL int
    fmt.Scanln(&numTestsL)
    testExtendedEuclideanLength(numTestsL)
}
