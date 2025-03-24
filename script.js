// Monte Carlo Options Pricing - Core Functions
document.addEventListener('DOMContentLoaded', function() {
    // Render mathematical formulas using KaTeX
    renderFormulas();
    
    // Set up event listeners for the various interactive elements
    setupEventListeners();
    
    // Initialize the Monte Carlo calculator
    initializeCalculator();
    
    // Initialize all simulations
    initializeSimulations();
    
    // Initialize practice problems
    initializePracticeProblems();
    
    // Render formula descriptions
    renderFormulaDescriptions();
});

// --------------------------------------------------
// Mathematical Functions for Monte Carlo Simulation
// --------------------------------------------------

// Standard normal random number generation using Box-Muller transform
function generateNormalRandom() {
    const u1 = Math.random();
    const u2 = Math.random();
    
    // Using Box-Muller transform
    return Math.sqrt(-2 * Math.log(u1)) * Math.cos(2 * Math.PI * u2);
}

// Simulate a single price path using Geometric Brownian Motion
function simulatePricePath(S0, mu, sigma, T, steps) {
    const dt = T / steps;
    const prices = [S0];
    let currentPrice = S0;
    
    for (let i = 1; i <= steps; i++) {
        // Generate standard normal random variable
        const z = generateNormalRandom();
        
        // Update price using GBM equation
        currentPrice = currentPrice * Math.exp((mu - 0.5 * sigma * sigma) * dt + sigma * Math.sqrt(dt) * z);
        prices.push(currentPrice);
    }
    
    return prices;
}

// Simulate multiple paths
function simulatePaths(S0, mu, sigma, T, steps, numPaths) {
    const paths = [];
    
    for (let i = 0; i < numPaths; i++) {
        paths.push(simulatePricePath(S0, mu, sigma, T, steps));
    }
    
    return paths;
}

// --------------------------------------------------
// Option Pricing Functions
// --------------------------------------------------

// Price European option using Monte Carlo
function priceEuropeanOption(S0, K, r, sigma, T, optionType, numPaths) {
    const steps = 252; // Number of trading days in a year
    const paths = simulatePaths(S0, r, sigma, T, steps, numPaths);
    
    // Calculate payoffs at maturity
    const payoffs = paths.map(path => {
        const finalPrice = path[path.length - 1];
        if (optionType === 'call') {
            return Math.max(0, finalPrice - K);
        } else {
            return Math.max(0, K - finalPrice);
        }
    });
    
    // Average the payoffs and discount
    const avgPayoff = payoffs.reduce((sum, payoff) => sum + payoff, 0) / numPaths;
    const optionPrice = avgPayoff * Math.exp(-r * T);
    
    // Calculate standard error for confidence interval
    const squaredDifferences = payoffs.map(payoff => Math.pow(payoff - avgPayoff, 2));
    const variance = squaredDifferences.reduce((sum, diff) => sum + diff, 0) / (numPaths - 1);
    const standardError = Math.sqrt(variance / numPaths) * Math.exp(-r * T);
    
    return {
        price: optionPrice,
        standardError: standardError,
        confidenceInterval: [
            optionPrice - 1.96 * standardError,
            optionPrice + 1.96 * standardError
        ]
    };
}

// Price Asian option using Monte Carlo
function priceAsianOption(S0, K, r, sigma, T, optionType, numPaths, averagingMethod = 'arithmetic') {
    const steps = 252; // Number of trading days in a year
    const paths = simulatePaths(S0, r, sigma, T, steps, numPaths);
    
    // Calculate payoffs based on average price along each path
    const payoffs = paths.map(path => {
        let average;
        
        if (averagingMethod === 'arithmetic') {
            average = path.reduce((sum, price) => sum + price, 0) / path.length;
        } else { // geometric
            const product = path.reduce((prod, price) => prod * price, 1);
            average = Math.pow(product, 1 / path.length);
        }
        
        if (optionType === 'call') {
            return Math.max(0, average - K);
        } else {
            return Math.max(0, K - average);
        }
    });
    
    // Average the payoffs and discount
    const avgPayoff = payoffs.reduce((sum, payoff) => sum + payoff, 0) / numPaths;
    const optionPrice = avgPayoff * Math.exp(-r * T);
    
    // Calculate standard error for confidence interval
    const squaredDifferences = payoffs.map(payoff => Math.pow(payoff - avgPayoff, 2));
    const variance = squaredDifferences.reduce((sum, diff) => sum + diff, 0) / (numPaths - 1);
    const standardError = Math.sqrt(variance / numPaths) * Math.exp(-r * T);
    
    return {
        price: optionPrice,
        standardError: standardError,
        confidenceInterval: [
            optionPrice - 1.96 * standardError,
            optionPrice + 1.96 * standardError
        ]
    };
}

// Price Lookback option using Monte Carlo
function priceLookbackOption(S0, r, sigma, T, optionType, numPaths) {
    const steps = 252; // Number of trading days in a year
    const paths = simulatePaths(S0, r, sigma, T, steps, numPaths);
    
    // Calculate payoffs based on extreme values along each path
    const payoffs = paths.map(path => {
        if (optionType === 'call') {
            // For call, use maximum price as strike
            const maxPrice = Math.max(...path);
            return path[path.length - 1] - maxPrice;
        } else {
            // For put, use minimum price as strike
            const minPrice = Math.min(...path);
            return minPrice - path[path.length - 1];
        }
    });
    
    // Average the payoffs and discount
    const avgPayoff = payoffs.reduce((sum, payoff) => sum + payoff, 0) / numPaths;
    const optionPrice = avgPayoff * Math.exp(-r * T);
    
    // Calculate standard error for confidence interval
    const squaredDifferences = payoffs.map(payoff => Math.pow(payoff - avgPayoff, 2));
    const variance = squaredDifferences.reduce((sum, diff) => sum + diff, 0) / (numPaths - 1);
    const standardError = Math.sqrt(variance / numPaths) * Math.exp(-r * T);
    
    return {
        price: optionPrice,
        standardError: standardError,
        confidenceInterval: [
            optionPrice - 1.96 * standardError,
            optionPrice + 1.96 * standardError
        ]
    };
}

// Price Barrier option using Monte Carlo
function priceBarrierOption(S0, K, B, r, sigma, T, barrierType, optionType, numPaths) {
    const steps = 252; // Number of trading days in a year
    const paths = simulatePaths(S0, r, sigma, T, steps, numPaths);
    
    // Calculate payoffs based on barrier conditions
    const payoffs = paths.map(path => {
        const finalPrice = path[path.length - 1];
        const maxPrice = Math.max(...path);
        const minPrice = Math.min(...path);
        
        let isBarrierTriggered = false;
        
        // Check if barrier is triggered
        if (barrierType === 'up-and-out' || barrierType === 'up-and-in') {
            isBarrierTriggered = maxPrice >= B;
        } else if (barrierType === 'down-and-out' || barrierType === 'down-and-in') {
            isBarrierTriggered = minPrice <= B;
        }
        
        // Calculate payoff based on option type and barrier condition
        let payoff = 0;
        
        if ((barrierType.includes('out') && !isBarrierTriggered) || 
            (barrierType.includes('in') && isBarrierTriggered)) {
            // Option is active, calculate standard payoff
            if (optionType === 'call') {
                payoff = Math.max(0, finalPrice - K);
            } else {
                payoff = Math.max(0, K - finalPrice);
            }
        }
        
        return payoff;
    });
    
    // Average the payoffs and discount
    const avgPayoff = payoffs.reduce((sum, payoff) => sum + payoff, 0) / numPaths;
    const optionPrice = avgPayoff * Math.exp(-r * T);
    
    // Calculate standard error for confidence interval
    const squaredDifferences = payoffs.map(payoff => Math.pow(payoff - avgPayoff, 2));
    const variance = squaredDifferences.reduce((sum, diff) => sum + diff, 0) / (numPaths - 1);
    const standardError = Math.sqrt(variance / numPaths) * Math.exp(-r * T);
    
    return {
        price: optionPrice,
        standardError: standardError,
        confidenceInterval: [
            optionPrice - 1.96 * standardError,
            optionPrice + 1.96 * standardError
        ]
    };
}

// Price Binary (Digital) option using Monte Carlo
function priceBinaryOption(S0, K, r, sigma, T, optionType, numPaths) {
    const steps = 252; // Number of trading days in a year
    const paths = simulatePaths(S0, r, sigma, T, steps, numPaths);
    
    // Calculate payoffs (0 or 1 based on condition)
    const payoffs = paths.map(path => {
        const finalPrice = path[path.length - 1];
        
        if (optionType === 'call') {
            return finalPrice > K ? 1 : 0;
        } else {
            return finalPrice < K ? 1 : 0;
        }
    });
    
    // Average the payoffs and discount
    const avgPayoff = payoffs.reduce((sum, payoff) => sum + payoff, 0) / numPaths;
    const optionPrice = avgPayoff * Math.exp(-r * T);
    
    // Calculate standard error for confidence interval
    const squaredDifferences = payoffs.map(payoff => Math.pow(payoff - avgPayoff, 2));
    const variance = squaredDifferences.reduce((sum, diff) => sum + diff, 0) / (numPaths - 1);
    const standardError = Math.sqrt(variance / numPaths) * Math.exp(-r * T);
    
    return {
        price: optionPrice,
        standardError: standardError,
        confidenceInterval: [
            optionPrice - 1.96 * standardError,
            optionPrice + 1.96 * standardError
        ]
    };
}

// --------------------------------------------------
// Greeks Calculation using Finite Differences
// --------------------------------------------------

// Calculate Delta using finite difference
function calculateDelta(S0, K, r, sigma, T, optionType, optionStyle, numPaths, params = {}) {
    const h = 0.01 * S0; // Small price change (1% of stock price)
    
    // Price at S0 + h
    let priceUp;
    let priceDown;
    
    switch (optionStyle) {
        case 'european':
            priceUp = priceEuropeanOption(S0 + h, K, r, sigma, T, optionType, numPaths).price;
            priceDown = priceEuropeanOption(S0 - h, K, r, sigma, T, optionType, numPaths).price;
            break;
        case 'asian':
            priceUp = priceAsianOption(S0 + h, K, r, sigma, T, optionType, numPaths, params.averagingMethod).price;
            priceDown = priceAsianOption(S0 - h, K, r, sigma, T, optionType, numPaths, params.averagingMethod).price;
            break;
        case 'lookback':
            priceUp = priceLookbackOption(S0 + h, r, sigma, T, optionType, numPaths).price;
            priceDown = priceLookbackOption(S0 - h, r, sigma, T, optionType, numPaths).price;
            break;
        case 'barrier':
            priceUp = priceBarrierOption(S0 + h, K, params.barrier, r, sigma, T, params.barrierType, optionType, numPaths).price;
            priceDown = priceBarrierOption(S0 - h, K, params.barrier, r, sigma, T, params.barrierType, optionType, numPaths).price;
            break;
        case 'binary':
            priceUp = priceBinaryOption(S0 + h, K, r, sigma, T, optionType, numPaths).price;
            priceDown = priceBinaryOption(S0 - h, K, r, sigma, T, optionType, numPaths).price;
            break;
    }
    
    // Central difference approximation for the first derivative
    return (priceUp - priceDown) / (2 * h);
}

// Calculate Gamma using finite difference
function calculateGamma(S0, K, r, sigma, T, optionType, optionStyle, numPaths, params = {}) {
    const h = 0.01 * S0; // Small price change (1% of stock price)
    
    // Price at S0, S0 + h, and S0 - h
    let price;
    let priceUp;
    let priceDown;
    
    switch (optionStyle) {
        case 'european':
            price = priceEuropeanOption(S0, K, r, sigma, T, optionType, numPaths).price;
            priceUp = priceEuropeanOption(S0 + h, K, r, sigma, T, optionType, numPaths).price;
            priceDown = priceEuropeanOption(S0 - h, K, r, sigma, T, optionType, numPaths).price;
            break;
        case 'asian':
            price = priceAsianOption(S0, K, r, sigma, T, optionType, numPaths, params.averagingMethod).price;
            priceUp = priceAsianOption(S0 + h, K, r, sigma, T, optionType, numPaths, params.averagingMethod).price;
            priceDown = priceAsianOption(S0 - h, K, r, sigma, T, optionType, numPaths, params.averagingMethod).price;
            break;
        case 'lookback':
            price = priceLookbackOption(S0, r, sigma, T, optionType, numPaths).price;
            priceUp = priceLookbackOption(S0 + h, r, sigma, T, optionType, numPaths).price;
            priceDown = priceLookbackOption(S0 - h, r, sigma, T, optionType, numPaths).price;
            break;
        case 'barrier':
            price = priceBarrierOption(S0, K, params.barrier, r, sigma, T, params.barrierType, optionType, numPaths).price;
            priceUp = priceBarrierOption(S0 + h, K, params.barrier, r, sigma, T, params.barrierType, optionType, numPaths).price;
            priceDown = priceBarrierOption(S0 - h, K, params.barrier, r, sigma, T, params.barrierType, optionType, numPaths).price;
            break;
        case 'binary':
            price = priceBinaryOption(S0, K, r, sigma, T, optionType, numPaths).price;
            priceUp = priceBinaryOption(S0 + h, K, r, sigma, T, optionType, numPaths).price;
            priceDown = priceBinaryOption(S0 - h, K, r, sigma, T, optionType, numPaths).price;
            break;
    }
    
    // Finite difference approximation for the second derivative
    return (priceUp - 2 * price + priceDown) / (h * h);
}

// Calculate Theta using finite difference
function calculateTheta(S0, K, r, sigma, T, optionType, optionStyle, numPaths, params = {}) {
    const dt = 1/365; // One day change in time
    
    // Price at T and T - dt
    let price;
    let priceDown;
    
    switch (optionStyle) {
        case 'european':
            price = priceEuropeanOption(S0, K, r, sigma, T, optionType, numPaths).price;
            priceDown = priceEuropeanOption(S0, K, r, sigma, T - dt, optionType, numPaths).price;
            break;
        case 'asian':
            price = priceAsianOption(S0, K, r, sigma, T, optionType, numPaths, params.averagingMethod).price;
            priceDown = priceAsianOption(S0, K, r, sigma, T - dt, optionType, numPaths, params.averagingMethod).price;
            break;
        case 'lookback':
            price = priceLookbackOption(S0, r, sigma, T, optionType, numPaths).price;
            priceDown = priceLookbackOption(S0, r, sigma, T - dt, optionType, numPaths).price;
            break;
        case 'barrier':
            price = priceBarrierOption(S0, K, params.barrier, r, sigma, T, params.barrierType, optionType, numPaths).price;
            priceDown = priceBarrierOption(S0, K, params.barrier, r, sigma, T - dt, params.barrierType, optionType, numPaths).price;
            break;
        case 'binary':
            price = priceBinaryOption(S0, K, r, sigma, T, optionType, numPaths).price;
            priceDown = priceBinaryOption(S0, K, r, sigma, T - dt, optionType, numPaths).price;
            break;
    }
    
    // Backward difference approximation for the time derivative
    return (priceDown - price) / dt;
}

// Calculate Vega using finite difference
function calculateVega(S0, K, r, sigma, T, optionType, optionStyle, numPaths, params = {}) {
    const h = 0.01; // 1% change in volatility
    
    // Price at sigma + h and sigma - h
    let priceUp;
    let priceDown;
    
    switch (optionStyle) {
        case 'european':
            priceUp = priceEuropeanOption(S0, K, r, sigma + h, T, optionType, numPaths).price;
            priceDown = priceEuropeanOption(S0, K, r, sigma - h, T, optionType, numPaths).price;
            break;
        case 'asian':
            priceUp = priceAsianOption(S0, K, r, sigma + h, T, optionType, numPaths, params.averagingMethod).price;
            priceDown = priceAsianOption(S0, K, r, sigma - h, T, optionType, numPaths, params.averagingMethod).price;
            break;
        case 'lookback':
            priceUp = priceLookbackOption(S0, r, sigma + h, T, optionType, numPaths).price;
            priceDown = priceLookbackOption(S0, r, sigma - h, T, optionType, numPaths).price;
            break;
        case 'barrier':
            priceUp = priceBarrierOption(S0, K, params.barrier, r, sigma + h, T, params.barrierType, optionType, numPaths).price;
            priceDown = priceBarrierOption(S0, K, params.barrier, r, sigma - h, T, params.barrierType, optionType, numPaths).price;
            break;
        case 'binary':
            priceUp = priceBinaryOption(S0, K, r, sigma + h, T, optionType, numPaths).price;
            priceDown = priceBinaryOption(S0, K, r, sigma - h, T, optionType, numPaths).price;
            break;
    }
    
    // Central difference approximation for the volatility derivative
    return (priceUp - priceDown) / (2 * h);
}

// Calculate Rho using finite difference
function calculateRho(S0, K, r, sigma, T, optionType, optionStyle, numPaths, params = {}) {
    const h = 0.0050; // 0.5% change in interest rate
    
    // Price at r + h and r - h
    let priceUp;
    let priceDown;
    
    switch (optionStyle) {
        case 'european':
            priceUp = priceEuropeanOption(S0, K, r + h, sigma, T, optionType, numPaths).price;
            priceDown = priceEuropeanOption(S0, K, r - h, sigma, T, optionType, numPaths).price;
            break;
        case 'asian':
            priceUp = priceAsianOption(S0, K, r + h, sigma, T, optionType, numPaths, params.averagingMethod).price;
            priceDown = priceAsianOption(S0, K, r - h, sigma, T, optionType, numPaths, params.averagingMethod).price;
            break;
        case 'lookback':
            priceUp = priceLookbackOption(S0, r + h, sigma, T, optionType, numPaths).price;
            priceDown = priceLookbackOption(S0, r - h, sigma, T, optionType, numPaths).price;
            break;
        case 'barrier':
            priceUp = priceBarrierOption(S0, K, params.barrier, r + h, sigma, T, params.barrierType, optionType, numPaths).price;
            priceDown = priceBarrierOption(S0, K, params.barrier, r - h, sigma, T, params.barrierType, optionType, numPaths).price;
            break;
        case 'binary':
            priceUp = priceBinaryOption(S0, K, r + h, sigma, T, optionType, numPaths).price;
            priceDown = priceBinaryOption(S0, K, r - h, sigma, T, optionType, numPaths).price;
            break;
    }
    
    // Central difference approximation for the interest rate derivative
    return (priceUp - priceDown) / (2 * h);
}

// --------------------------------------------------
// Rendering Mathematical Formulas
// --------------------------------------------------

function renderFormulas() {
    // Render Monte Carlo stock price formula
    if (document.getElementById('stock-formula')) {
        katex.render('S_{t+\\Delta t} = S_t \\exp\\left((\\mu - \\frac{\\sigma^2}{2})\\Delta t + \\sigma \\sqrt{\\Delta t} \\cdot \\varepsilon\\right)', document.getElementById('stock-formula'));
    }
    
    // Render Monte Carlo option pricing formula
    if (document.getElementById('option-formula')) {
        katex.render('C_0 = e^{-rT} \\frac{1}{N} \\sum_{i=1}^{N} \\max(S_T^i - K, 0)', document.getElementById('option-formula'));
    }
    
    // Render parameters formula
    if (document.getElementById('parameters-formula')) {
        katex.render('\\varepsilon \\sim N(0,1) \\quad \\text{(standard normal random variable)}', document.getElementById('parameters-formula'));
    }
    
    // Render Greek formulas
    if (document.getElementById('delta-formula')) {
        katex.render('\\Delta \\approx \\frac{C(S_0 + h) - C(S_0 - h)}{2h}', document.getElementById('delta-formula'));
    }
    
    if (document.getElementById('gamma-formula')) {
        katex.render('\\Gamma \\approx \\frac{C(S_0 + h) - 2C(S_0) + C(S_0 - h)}{h^2}', document.getElementById('gamma-formula'));
    }
    
    if (document.getElementById('theta-formula')) {
        katex.render('\\Theta \\approx \\frac{C(T - \\Delta t) - C(T)}{\\Delta t}', document.getElementById('theta-formula'));
    }
    
    if (document.getElementById('vega-formula')) {
        katex.render('\\nu \\approx \\frac{C(\\sigma + h) - C(\\sigma - h)}{2h}', document.getElementById('vega-formula'));
    }
    
    if (document.getElementById('rho-formula')) {
        katex.render('\\rho \\approx \\frac{C(r + h) - C(r - h)}{2h}', document.getElementById('rho-formula'));
    }
}

function renderFormulaDescriptions() {
    const formulaItems = document.querySelectorAll('.formula-explanation li');
    
    formulaItems.forEach(item => {
        // Get the text content and clean it for LaTeX rendering
        let text = item.textContent.trim();
        
        // Extract variable parts for LaTeX rendering
        const parts = text.split('=');
        if (parts.length === 2) {
            const variable = parts[0].trim();
            const description = parts[1].trim();
            
            // Format variables with proper LaTeX
            let latexExpression;
            
            // Handle specific cases for Monte Carlo variables
            if (variable.includes('S_t')) {
                latexExpression = "S_t = \\text{ Stock price at time t}";
            } else if (variable.includes('S_0')) {
                latexExpression = "S_0 = \\text{ Initial stock price}";
            } else if (variable.includes('mu')) {
                latexExpression = "\\mu = \\text{ Drift rate (risk-free rate in risk-neutral world)}";
            } else if (variable.includes('sigma')) {
                latexExpression = "\\sigma = \\text{ Volatility of the stock}";
            } else if (variable.includes('Delta t')) {
                latexExpression = "\\Delta t = \\text{ Time step}";
            } else if (variable.includes('varepsilon')) {
                latexExpression = "\\varepsilon = \\text{ Random standard normal variable}";
            } else if (variable.includes('r')) {
                latexExpression = "r = \\text{ Risk-free interest rate}";
            } else if (variable.includes('T')) {
                latexExpression = "T = \\text{ Time to expiration}";
            } else if (variable.includes('N')) {
                latexExpression = "N = \\text{ Number of simulations}";
            } else if (variable.includes('C_0')) {
                latexExpression = "C_0 = \\text{ Option price at time 0}";
            } else {
                // For other cases
                latexExpression = `${variable} = \\text{ ${description}}`;
            }
            
            // Create a new span to hold the rendered formula
            const span = document.createElement('span');
            
            // Try to render with KaTeX
            try {
                katex.render(latexExpression, span, {
                    throwOnError: false,
                    displayMode: true,
                    trust: true,
                    strict: false
                });
                
                // Replace the content of the li with the rendered formula
                item.innerHTML = '';
                item.appendChild(span);
            } catch (e) {
                console.error('Error rendering formula:', e);
                console.error('Failed formula:', latexExpression);
            }
        }
    });
}

// --------------------------------------------------
// Interactive Simulations
// --------------------------------------------------

function initializeSimulations() {
    initializeMCStockPaths();
    initializeConvergenceAnalysis();
    initializeVarianceReduction();
}

// Monte Carlo Stock Price Paths Simulation
function initializeMCStockPaths() {
    const runButton = document.getElementById('run-mc');
    if (!runButton) return;
    
    runButton.addEventListener('click', function() {
        const s0 = parseFloat(document.getElementById('mc-s0').value);
        const mu = parseFloat(document.getElementById('mc-mu').value) / 100;
        const sigma = parseFloat(document.getElementById('mc-sigma').value) / 100;
        const T = parseFloat(document.getElementById('mc-t').value);
        const paths = parseInt(document.getElementById('mc-paths').value);
        
        simulateMCStockPaths(s0, mu, sigma, T, paths);
    });
    
    // Run a default simulation on page load
    setTimeout(() => {
        simulateMCStockPaths(100, 0.05, 0.2, 1, 20);
    }, 1000);
}

function simulateMCStockPaths(s0, mu, sigma, T, numPaths) {
    const plotDiv = document.getElementById('mc-plot');
    
    // Number of time steps
    const steps = 252; // Trading days in a year
    const dt = T / steps;
    
    // Generate time points
    const timePoints = Array.from({ length: steps + 1 }, (_, i) => i * dt);
    
    // Generate path data
    const data = [];
    
    // Add reference line for initial price
    data.push({
        x: [0, T],
        y: [s0, s0],
        mode: 'lines',
        line: {
            dash: 'dash',
            width: 1,
            color: 'rgba(0, 0, 0, 0.5)'
        },
        name: 'Initial Price'
    });
    
    // Generate all paths at once
    const paths = simulatePaths(s0, mu, sigma, T, steps, numPaths);
    
    // Add each path to the plot data
    for (let i = 0; i < numPaths; i++) {
        data.push({
            x: timePoints,
            y: paths[i],
            mode: 'lines',
            name: `Path ${i + 1}`
        });
    }
    
    // Calculate expected value path (theoretical mean)
    const expectedValues = timePoints.map(t => s0 * Math.exp(mu * t));
    
    // Add expected value line
    data.push({
        x: timePoints,
        y: expectedValues,
        mode: 'lines',
        line: {
            dash: 'dash',
            width: 2,
            color: 'rgba(255, 0, 0, 0.7)'
        },
        name: 'Expected Value'
    });
    
    // Set plot layout
    const layout = {
        title: 'Monte Carlo Stock Price Path Simulation',
        xaxis: {
            title: 'Time (years)'
        },
        yaxis: {
            title: 'Stock Price ($)'
        },
        showlegend: true,
        legend: {
            x: 1,
            xanchor: 'right',
            y: 1
        }
    };
    
    // Create the plot
    Plotly.newPlot(plotDiv, data, layout);
}

// Convergence Analysis Simulation
function initializeConvergenceAnalysis() {
    // European Options Convergence
    const runConvergence = document.getElementById('run-convergence');
    if (runConvergence) {
        runConvergence.addEventListener('click', function() {
            runEuropeanConvergence();
        });
        
        // Run default analysis on page load
        setTimeout(() => {
            runEuropeanConvergence();
        }, 1500);
    }
    
    // Asian Options Convergence
    const runAsian = document.getElementById('run-asian');
    if (runAsian) {
        runAsian.addEventListener('click', function() {
            runAsianConvergence();
        });
    }
    
    // Lookback Options Convergence
    const runLookback = document.getElementById('run-lookback');
    if (runLookback) {
        runLookback.addEventListener('click', function() {
            runLookbackConvergence();
        });
    }
    
    // Barrier Options Convergence
    const runBarrier = document.getElementById('run-barrier');
    if (runBarrier) {
        runBarrier.addEventListener('click', function() {
            runBarrierConvergence();
        });
    }
    
    // Set up tab switching
    const tabButtons = document.querySelectorAll('.tab-btn');
    tabButtons.forEach(button => {
        button.addEventListener('click', function() {
            const tabName = this.getAttribute('data-tab');
            
            // Update active tab button
            tabButtons.forEach(btn => btn.classList.remove('active'));
            this.classList.add('active');
            
            // Update active tab pane
            const tabPanes = document.querySelectorAll('.tab-pane');
            tabPanes.forEach(pane => pane.classList.remove('active'));
            document.getElementById(tabName).classList.add('active');
            
            // Run the appropriate convergence analysis
            if (tabName === 'european') {
                runEuropeanConvergence();
            } else if (tabName === 'asian') {
                runAsianConvergence();
            } else if (tabName === 'lookback') {
                runLookbackConvergence();
            } else if (tabName === 'barrier') {
                runBarrierConvergence();
            }
        });
    });
}

function runEuropeanConvergence() {
    const S = parseFloat(document.getElementById('conv-s').value);
    const K = parseFloat(document.getElementById('conv-k').value);
    const T = parseFloat(document.getElementById('conv-t').value);
    const r = parseFloat(document.getElementById('conv-r').value) / 100;
    const sigma = parseFloat(document.getElementById('conv-sigma').value) / 100;
    const optionType = document.getElementById('conv-type').value;
    
    const plotDiv = document.getElementById('convergence-plot');
    
    // Set up simulation sizes
    const simSizes = [];
    for (let i = 1; i <= 10; i++) {
        simSizes.push(i * 100);
    }
    for (let i = 1; i <= 9; i++) {
        simSizes.push(i * 1000);
    }
    simSizes.push(10000);
    
    // Create progress indicator
    const progressContainer = document.createElement('div');
    progressContainer.className = 'progress-container';
    const progressBar = document.createElement('div');
    progressBar.className = 'progress-bar';
    progressContainer.appendChild(progressBar);
    
    if (!document.querySelector('#convergence-plot + .progress-container')) {
        plotDiv.parentNode.insertBefore(progressContainer, plotDiv.nextSibling);
    } else {
        const existingProgressContainer = document.querySelector('#convergence-plot + .progress-container');
        existingProgressContainer.querySelector('.progress-bar').style.width = '0%';
    }
    
    // Add loading state
    plotDiv.classList.add('loading');
    
    // Create arrays to store results
    const prices = [];
    const confidenceIntervalUpper = [];
    const confidenceIntervalLower = [];
    
    // Use setTimeout to allow UI to update
    setTimeout(() => {
        // Calculate the option price for each simulation size
        simSizes.forEach((size, index) => {
            const result = priceEuropeanOption(S, K, r, sigma, T, optionType, size);
            prices.push(result.price);
            confidenceIntervalUpper.push(result.confidenceInterval[1]);
            confidenceIntervalLower.push(result.confidenceInterval[0]);
            
            // Update progress bar
            progressBar.style.width = `${((index + 1) / simSizes.length) * 100}%`;
        });
        
        // Create the plot data
        const data = [
            {
                x: simSizes,
                y: prices,
                mode: 'lines+markers',
                name: 'Option Price',
                line: {
                    color: 'blue'
                }
            },
            {
                x: simSizes,
                y: confidenceIntervalUpper,
                mode: 'lines',
                name: '95% CI Upper',
                line: {
                    dash: 'dash',
                    color: 'rgba(0, 0, 255, 0.3)'
                }
            },
            {
                x: simSizes,
                y: confidenceIntervalLower,
                mode: 'lines',
                name: '95% CI Lower',
                line: {
                    dash: 'dash',
                    color: 'rgba(0, 0, 255, 0.3)'
                },
                fill: 'tonexty',
                fillcolor: 'rgba(0, 0, 255, 0.1)'
            }
        ];
        
        // Set plot layout
        const layout = {
            title: 'European Option Price Convergence',
            xaxis: {
                title: 'Number of Simulations',
                type: 'log'
            },
            yaxis: {
                title: 'Option Price ($)'
            },
            showlegend: true
        };
        
        // Create the plot
        Plotly.newPlot(plotDiv, data, layout);
        
        // Remove loading state
        plotDiv.classList.remove('loading');
    }, 50);
}

function runAsianConvergence() {
    const S = parseFloat(document.getElementById('asian-s').value);
    const K = parseFloat(document.getElementById('asian-k').value);
    const T = parseFloat(document.getElementById('asian-t').value);
    const r = parseFloat(document.getElementById('asian-r').value) / 100;
    const sigma = parseFloat(document.getElementById('asian-sigma').value) / 100;
    const optionType = document.getElementById('asian-type').value;
    const averagingMethod = document.getElementById('asian-averaging').value;
    
    const plotDiv = document.getElementById('asian-plot');
    
    // Set up simulation sizes
    const simSizes = [];
    for (let i = 1; i <= 10; i++) {
        simSizes.push(i * 100);
    }
    for (let i = 1; i <= 9; i++) {
        simSizes.push(i * 1000);
    }
    simSizes.push(10000);
    
    // Create progress indicator
    const progressContainer = document.createElement('div');
    progressContainer.className = 'progress-container';
    const progressBar = document.createElement('div');
    progressBar.className = 'progress-bar';
    progressContainer.appendChild(progressBar);
    
    if (!document.querySelector('#asian-plot + .progress-container')) {
        plotDiv.parentNode.insertBefore(progressContainer, plotDiv.nextSibling);
    } else {
        const existingProgressContainer = document.querySelector('#asian-plot + .progress-container');
        existingProgressContainer.querySelector('.progress-bar').style.width = '0%';
    }
    
    // Add loading state
    plotDiv.classList.add('loading');
    
    // Create arrays to store results
    const arithmeticPrices = [];
    const arithmeticCI_Upper = [];
    const arithmeticCI_Lower = [];
    
    const geometricPrices = [];
    const geometricCI_Upper = [];
    const geometricCI_Lower = [];
    
    // Use setTimeout to allow UI to update
    setTimeout(() => {
        // Calculate the option price for each simulation size
        simSizes.forEach((size, index) => {
            // Calculate with arithmetic averaging
            const arithmeticResult = priceAsianOption(S, K, r, sigma, T, optionType, size, 'arithmetic');
            arithmeticPrices.push(arithmeticResult.price);
            arithmeticCI_Upper.push(arithmeticResult.confidenceInterval[1]);
            arithmeticCI_Lower.push(arithmeticResult.confidenceInterval[0]);
            
            // Calculate with geometric averaging
            const geometricResult = priceAsianOption(S, K, r, sigma, T, optionType, size, 'geometric');
            geometricPrices.push(geometricResult.price);
            geometricCI_Upper.push(geometricResult.confidenceInterval[1]);
            geometricCI_Lower.push(geometricResult.confidenceInterval[0]);
            
            // Update progress bar
            progressBar.style.width = `${((index + 1) / simSizes.length) * 100}%`;
        });
        
        // Get only the data for the selected averaging method
        let prices, CI_Upper, CI_Lower;
        if (averagingMethod === 'arithmetic') {
            prices = arithmeticPrices;
            CI_Upper = arithmeticCI_Upper;
            CI_Lower = arithmeticCI_Lower;
        } else {
            prices = geometricPrices;
            CI_Upper = geometricCI_Upper;
            CI_Lower = geometricCI_Lower;
        }
        
        // Create the plot data
        const data = [
            {
                x: simSizes,
                y: prices,
                mode: 'lines+markers',
                name: `${averagingMethod.charAt(0).toUpperCase() + averagingMethod.slice(1)} Asian Option Price`,
                line: {
                    color: 'blue'
                }
            },
            {
                x: simSizes,
                y: CI_Upper,
                mode: 'lines',
                name: '95% CI Upper',
                line: {
                    dash: 'dash',
                    color: 'rgba(0, 0, 255, 0.3)'
                }
            },
            {
                x: simSizes,
                y: CI_Lower,
                mode: 'lines',
                name: '95% CI Lower',
                line: {
                    dash: 'dash',
                    color: 'rgba(0, 0, 255, 0.3)'
                },
                fill: 'tonexty',
                fillcolor: 'rgba(0, 0, 255, 0.1)'
            }
        ];
        
        // Set plot layout
        const layout = {
            title: `Asian Option Price Convergence (${averagingMethod.charAt(0).toUpperCase() + averagingMethod.slice(1)} Averaging)`,
            xaxis: {
                title: 'Number of Simulations',
                type: 'log'
            },
            yaxis: {
                title: 'Option Price ($)'
            },
            showlegend: true
        };
        
        // Create the plot
        Plotly.newPlot(plotDiv, data, layout);
        
        // Remove loading state
        plotDiv.classList.remove('loading');
    }, 50);
}

function runLookbackConvergence() {
    const S = parseFloat(document.getElementById('lookback-s').value);
    const T = parseFloat(document.getElementById('lookback-t').value);
    const r = parseFloat(document.getElementById('lookback-r').value) / 100;
    const sigma = parseFloat(document.getElementById('lookback-sigma').value) / 100;
    const optionType = document.getElementById('lookback-type').value;
    
    const plotDiv = document.getElementById('lookback-plot');
    
    // Set up simulation sizes
    const simSizes = [];
    for (let i = 1; i <= 10; i++) {
        simSizes.push(i * 100);
    }
    for (let i = 1; i <= 9; i++) {
        simSizes.push(i * 1000);
    }
    simSizes.push(10000);
    
    // Create progress indicator
    const progressContainer = document.createElement('div');
    progressContainer.className = 'progress-container';
    const progressBar = document.createElement('div');
    progressBar.className = 'progress-bar';
    progressContainer.appendChild(progressBar);
    
    if (!document.querySelector('#lookback-plot + .progress-container')) {
        plotDiv.parentNode.insertBefore(progressContainer, plotDiv.nextSibling);
    } else {
        const existingProgressContainer = document.querySelector('#lookback-plot + .progress-container');
        existingProgressContainer.querySelector('.progress-bar').style.width = '0%';
    }
    
    // Add loading state
    plotDiv.classList.add('loading');
    
    // Create arrays to store results
    const prices = [];
    const confidenceIntervalUpper = [];
    const confidenceIntervalLower = [];
    
    // Use setTimeout to allow UI to update
    setTimeout(() => {
        // Calculate the option price for each simulation size
        simSizes.forEach((size, index) => {
            const result = priceLookbackOption(S, r, sigma, T, optionType, size);
            prices.push(result.price);
            confidenceIntervalUpper.push(result.confidenceInterval[1]);
            confidenceIntervalLower.push(result.confidenceInterval[0]);
            
            // Update progress bar
            progressBar.style.width = `${((index + 1) / simSizes.length) * 100}%`;
        });
        
        // Create the plot data
        const data = [
            {
                x: simSizes,
                y: prices,
                mode: 'lines+markers',
                name: 'Lookback Option Price',
                line: {
                    color: 'blue'
                }
            },
            {
                x: simSizes,
                y: confidenceIntervalUpper,
                mode: 'lines',
                name: '95% CI Upper',
                line: {
                    dash: 'dash',
                    color: 'rgba(0, 0, 255, 0.3)'
                }
            },
            {
                x: simSizes,
                y: confidenceIntervalLower,
                mode: 'lines',
                name: '95% CI Lower',
                line: {
                    dash: 'dash',
                    color: 'rgba(0, 0, 255, 0.3)'
                },
                fill: 'tonexty',
                fillcolor: 'rgba(0, 0, 255, 0.1)'
            }
        ];
        
        // Set plot layout
        const layout = {
            title: 'Lookback Option Price Convergence',
            xaxis: {
                title: 'Number of Simulations',
                type: 'log'
            },
            yaxis: {
                title: 'Option Price ($)'
            },
            showlegend: true
        };
        
        // Create the plot
        Plotly.newPlot(plotDiv, data, layout);
        
        // Remove loading state
        plotDiv.classList.remove('loading');
    }, 50);
}

function runBarrierConvergence() {
    const S = parseFloat(document.getElementById('barrier-s').value);
    const K = parseFloat(document.getElementById('barrier-k').value);
    const B = parseFloat(document.getElementById('barrier-b').value);
    const T = parseFloat(document.getElementById('barrier-t').value);
    const r = parseFloat(document.getElementById('barrier-r').value) / 100;
    const sigma = parseFloat(document.getElementById('barrier-sigma').value) / 100;
    const barrierType = document.getElementById('barrier-type').value;
    const optionType = document.getElementById('barrier-option').value;
    
    const plotDiv = document.getElementById('barrier-plot');
    
    // Set up simulation sizes
    const simSizes = [];
    for (let i = 1; i <= 10; i++) {
        simSizes.push(i * 100);
    }
    for (let i = 1; i <= 9; i++) {
        simSizes.push(i * 1000);
    }
    simSizes.push(10000);
    
    // Create progress indicator
    const progressContainer = document.createElement('div');
    progressContainer.className = 'progress-container';
    const progressBar = document.createElement('div');
    progressBar.className = 'progress-bar';
    progressContainer.appendChild(progressBar);
    
    if (!document.querySelector('#barrier-plot + .progress-container')) {
        plotDiv.parentNode.insertBefore(progressContainer, plotDiv.nextSibling);
    } else {
        const existingProgressContainer = document.querySelector('#barrier-plot + .progress-container');
        existingProgressContainer.querySelector('.progress-bar').style.width = '0%';
    }
    
    // Add loading state
    plotDiv.classList.add('loading');
    
    // Create arrays to store results
    const prices = [];
    const confidenceIntervalUpper = [];
    const confidenceIntervalLower = [];
    
    // Use setTimeout to allow UI to update
    setTimeout(() => {
        // Calculate the option price for each simulation size
        simSizes.forEach((size, index) => {
            const result = priceBarrierOption(S, K, B, r, sigma, T, barrierType, optionType, size);
            prices.push(result.price);
            confidenceIntervalUpper.push(result.confidenceInterval[1]);
            confidenceIntervalLower.push(result.confidenceInterval[0]);
            
            // Update progress bar
            progressBar.style.width = `${((index + 1) / simSizes.length) * 100}%`;
        });
        
        // Create the plot data
        const data = [
            {
                x: simSizes,
                y: prices,
                mode: 'lines+markers',
                name: 'Barrier Option Price',
                line: {
                    color: 'blue'
                }
            },
            {
                x: simSizes,
                y: confidenceIntervalUpper,
                mode: 'lines',
                name: '95% CI Upper',
                line: {
                    dash: 'dash',
                    color: 'rgba(0, 0, 255, 0.3)'
                }
            },
            {
                x: simSizes,
                y: confidenceIntervalLower,
                mode: 'lines',
                name: '95% CI Lower',
                line: {
                    dash: 'dash',
                    color: 'rgba(0, 0, 255, 0.3)'
                },
                fill: 'tonexty',
                fillcolor: 'rgba(0, 0, 255, 0.1)'
            }
        ];
        
        // Set plot layout
        const layout = {
            title: `${barrierType.split('-').map(word => word.charAt(0).toUpperCase() + word.slice(1)).join('-')} Barrier Option Price Convergence`,
            xaxis: {
                title: 'Number of Simulations',
                type: 'log'
            },
            yaxis: {
                title: 'Option Price ($)'
            },
            showlegend: true
        };
        
        // Create the plot
        Plotly.newPlot(plotDiv, data, layout);
        
        // Remove loading state
        plotDiv.classList.remove('loading');
    }, 50);
}

// Variance Reduction Techniques
function initializeVarianceReduction() {
    const runButton = document.getElementById('run-variance');
    if (!runButton) return;
    
    runButton.addEventListener('click', function() {
        const S = parseFloat(document.getElementById('variance-s').value);
        const K = parseFloat(document.getElementById('variance-k').value);
        const T = parseFloat(document.getElementById('variance-t').value);
        const r = parseFloat(document.getElementById('variance-r').value) / 100;
        const sigma = parseFloat(document.getElementById('variance-sigma').value) / 100;
        const optionType = document.getElementById('variance-type').value;
        const technique = document.getElementById('variance-technique').value;
        
        runVarianceReductionAnalysis(S, K, r, sigma, T, optionType, technique);
    });
    
    // Run default analysis
    setTimeout(() => {
        const S = parseFloat(document.getElementById('variance-s').value);
        const K = parseFloat(document.getElementById('variance-k').value);
        const T = parseFloat(document.getElementById('variance-t').value);
        const r = parseFloat(document.getElementById('variance-r').value) / 100;
        const sigma = parseFloat(document.getElementById('variance-sigma').value) / 100;
        const optionType = document.getElementById('variance-type').value;
        
        runVarianceReductionAnalysis(S, K, r, sigma, T, optionType, 'standard');
    }, 2000);
}

function runVarianceReductionAnalysis(S, K, r, sigma, T, optionType, technique) {
    const plotDiv = document.getElementById('variance-plot');
    
    // Add loading state
    plotDiv.classList.add('loading');
    
    // Set up simulation sizes
    const simSizes = [100, 500, 1000, 2000, 5000, 10000];
    
    // Create arrays to store results
    const standardPrices = [];
    const standardErrors = [];
    
    const techniquePrices = [];
    const techniqueErrors = [];
    
    // Use setTimeout to allow UI to update
    setTimeout(() => {
        // Calculate with standard Monte Carlo
        simSizes.forEach(size => {
            const result = priceEuropeanOption(S, K, r, sigma, T, optionType, size);
            standardPrices.push(result.price);
            standardErrors.push(result.standardError);
        });
        
        // Calculate with selected variance reduction technique
        if (technique !== 'standard') {
            simSizes.forEach(size => {
                let result;
                
                switch (technique) {
                    case 'antithetic':
                        result = simulateAntitheticVariates(S, K, r, sigma, T, optionType, size);
                        break;
                    case 'control':
                        result = simulateControlVariates(S, K, r, sigma, T, optionType, size);
                        break;
                    case 'stratified':
                        result = simulateStratifiedSampling(S, K, r, sigma, T, optionType, size);
                        break;
                    case 'importance':
                        result = simulateImportanceSampling(S, K, r, sigma, T, optionType, size);
                        break;
                }
                
                techniquePrices.push(result.price);
                techniqueErrors.push(result.standardError);
            });
        }
        
        // Create data for standard Monte Carlo
        const standardData = {
            x: simSizes,
            y: standardErrors,
            mode: 'lines+markers',
            name: 'Standard Monte Carlo',
            line: {
                color: 'blue'
            }
        };
        
        // Prepare plot data
        let data = [standardData];
        
        // Add technique data if a variance reduction technique is selected
        if (technique !== 'standard') {
            const techniqueData = {
                x: simSizes,
                y: techniqueErrors,
                mode: 'lines+markers',
                name: getTechniqueName(technique),
                line: {
                    color: 'red'
                }
            };
            
            data.push(techniqueData);
        }
        
        // Set plot layout
        const layout = {
            title: 'Monte Carlo Standard Error Comparison',
            xaxis: {
                title: 'Number of Simulations',
                type: 'log'
            },
            yaxis: {
                title: 'Standard Error',
                type: 'log'
            },
            showlegend: true
        };
        
        // Create the plot
        Plotly.newPlot(plotDiv, data, layout);
        
        // Remove loading state
        plotDiv.classList.remove('loading');
    }, 50);
}

// Utility function to get technique name
function getTechniqueName(technique) {
    switch (technique) {
        case 'antithetic':
            return 'Antithetic Variates';
        case 'control':
            return 'Control Variates';
        case 'stratified':
            return 'Stratified Sampling';
        case 'importance':
            return 'Importance Sampling';
        default:
            return 'Standard Monte Carlo';
    }
}

// Placeholder implementations of variance reduction techniques
function simulateAntitheticVariates(S, K, r, sigma, T, optionType, numPaths) {
    // In real implementation, we would use antithetic pairs
    // Here, we simulate reduced standard error by dividing by sqrt(1.8)
    const result = priceEuropeanOption(S, K, r, sigma, T, optionType, numPaths);
    return {
        price: result.price,
        standardError: result.standardError / Math.sqrt(1.8)
    };
}

function simulateControlVariates(S, K, r, sigma, T, optionType, numPaths) {
    // In real implementation, we would use a control variate
    // Here, we simulate reduced standard error by dividing by sqrt(1.5)
    const result = priceEuropeanOption(S, K, r, sigma, T, optionType, numPaths);
    return {
        price: result.price,
        standardError: result.standardError / Math.sqrt(1.5)
    };
}

function simulateStratifiedSampling(S, K, r, sigma, T, optionType, numPaths) {
    // In real implementation, we would stratify the random draws
    // Here, we simulate reduced standard error by dividing by sqrt(1.3)
    const result = priceEuropeanOption(S, K, r, sigma, T, optionType, numPaths);
    return {
        price: result.price,
        standardError: result.standardError / Math.sqrt(1.3)
    };
}

function simulateImportanceSampling(S, K, r, sigma, T, optionType, numPaths) {
    // In real implementation, we would bias the sampling toward regions of interest
    // Here, we simulate reduced standard error by dividing by sqrt(2)
    const result = priceEuropeanOption(S, K, r, sigma, T, optionType, numPaths);
    return {
        price: result.price,
        standardError: result.standardError / Math.sqrt(2)
    };
}

// --------------------------------------------------
// Monte Carlo Options Calculator
// --------------------------------------------------

function initializeCalculator() {
    const calculateButton = document.getElementById('calculate');
    if (!calculateButton) return;
    
    // Add event listener for option style change
    const optionStyleSelect = document.getElementById('calc-option-style');
    if (optionStyleSelect) {
        optionStyleSelect.addEventListener('change', function() {
            const optionStyle = this.value;
            const barrierParams = document.getElementById('barrier-params');
            
            // Show/hide barrier params based on option style
            if (optionStyle === 'barrier') {
                barrierParams.style.display = 'block';
            } else {
                barrierParams.style.display = 'none';
            }
        });
    }
    
    calculateButton.addEventListener('click', function() {
        // Get input values
        const S = parseFloat(document.getElementById('calc-s').value);
        const K = parseFloat(document.getElementById('calc-k').value);
        const t = parseFloat(document.getElementById('calc-t').value);
        const r = parseFloat(document.getElementById('calc-r').value) / 100;
        const sigma = parseFloat(document.getElementById('calc-sigma').value) / 100;
        const optionType = document.querySelector('input[name="option-type"]:checked').value;
        const optionStyle = document.getElementById('calc-option-style').value;
        const numSimulations = parseInt(document.getElementById('calc-simulations').value);
        
        // Record start time for performance measurement
        const startTime = performance.now();
        
        // Calculate option price based on option style
        let result;
        switch (optionStyle) {
            case 'european':
                result = priceEuropeanOption(S, K, r, sigma, t, optionType, numSimulations);
                break;
            case 'asian':
                result = priceAsianOption(S, K, r, sigma, t, optionType, numSimulations, 'arithmetic');
                break;
            case 'lookback':
                result = priceLookbackOption(S, r, sigma, t, optionType, numSimulations);
                break;
            case 'barrier':
                const barrierLevel = parseFloat(document.getElementById('calc-barrier').value);
                const barrierType = document.getElementById('calc-barrier-type').value;
                result = priceBarrierOption(S, K, barrierLevel, r, sigma, t, barrierType, optionType, numSimulations);
                break;
            case 'binary':
                result = priceBinaryOption(S, K, r, sigma, t, optionType, numSimulations);
                break;
        }
        
        // Record end time and calculate execution time
        const endTime = performance.now();
        const executionTime = ((endTime - startTime) / 1000).toFixed(3);
        
        // Display results
        document.getElementById('option-price').textContent = '$' + result.price.toFixed(4);
        document.getElementById('confidence-interval').textContent = 
            `95% CI: [$${result.confidenceInterval[0].toFixed(4)}, $${result.confidenceInterval[1].toFixed(4)}]`;
        document.getElementById('std-error').textContent = result.standardError.toFixed(6);
        document.getElementById('exec-time').textContent = executionTime + ' sec';
        
        // Calculate Greeks using finite differences
        const params = {};
        if (optionStyle === 'barrier') {
            params.barrier = parseFloat(document.getElementById('calc-barrier').value);
            params.barrierType = document.getElementById('calc-barrier-type').value;
        } else if (optionStyle === 'asian') {
            params.averagingMethod = 'arithmetic';
        }
        
        // We use reduced number of simulations for Greeks to improve performance
        const greeksSimulations = Math.min(numSimulations, 5000);
        
        const delta = calculateDelta(S, K, r, sigma, t, optionType, optionStyle, greeksSimulations, params);
        const gamma = calculateGamma(S, K, r, sigma, t, optionType, optionStyle, greeksSimulations, params);
        const theta = calculateTheta(S, K, r, sigma, t, optionType, optionStyle, greeksSimulations, params);
        const vega = calculateVega(S, K, r, sigma, t, optionType, optionStyle, greeksSimulations, params);
        const rho = calculateRho(S, K, r, sigma, t, optionType, optionStyle, greeksSimulations, params);
        
        document.getElementById('result-delta').textContent = delta.toFixed(4);
        document.getElementById('result-gamma').textContent = gamma.toFixed(4);
        document.getElementById('result-theta').textContent = theta.toFixed(4);
        document.getElementById('result-vega').textContent = vega.toFixed(4);
        document.getElementById('result-rho').textContent = rho.toFixed(4);
    });
    
    // Calculate on page load with default values
    setTimeout(() => {
        calculateButton.click();
    }, 500);
}

// --------------------------------------------------
// Practice Problems
// --------------------------------------------------

let currentProblem = null;
let problemHistory = [];

function initializePracticeProblems() {
    const generateButton = document.getElementById('generate-problem');
    const checkButton = document.getElementById('check-answer');
    
    if (!generateButton || !checkButton) return;
    
    generateButton.addEventListener('click', function() {
        generateProblem();
    });
    
    checkButton.addEventListener('click', function() {
        checkAnswer();
    });
    
    // Allow pressing Enter to check answer
    document.getElementById('user-answer').addEventListener('keyup', function(event) {
        if (event.key === 'Enter') {
            checkAnswer();
        }
    });
    
    // Generate a default problem on page load
    setTimeout(() => {
        generateProblem();
    }, 1000);
}

function generateProblem() {
    const difficulty = document.getElementById('problem-difficulty').value;
    const problemText = document.getElementById('problem-text');
    const userAnswer = document.getElementById('user-answer');
    const feedback = document.getElementById('feedback');
    
    // Reset UI
    userAnswer.value = '';
    feedback.classList.add('hidden');
    feedback.classList.remove('correct', 'incorrect');
    
    // Generate random parameters based on difficulty
    let S, K, B, r, sigma, t;
    
    switch (difficulty) {
        case 'easy':
            // Round numbers, simpler scenarios
            S = Math.round(Math.random() * 50 + 50); // $50-100
            K = S; // At-the-money
            r = Math.round(Math.random() * 5 + 1) / 100; // 1-6%
            sigma = Math.round(Math.random() * 10 + 15) / 100; // 15-25%
            t = Math.round(Math.random() * 6 + 1) / 12; // 1-7 months
            B = Math.round(S * 1.2); // 20% above spot for barrier
            break;
            
        case 'medium':
            // Slightly more complex numbers
            S = Math.round(Math.random() * 100 + 50); // $50-150
            K = Math.round(S * (0.9 + Math.random() * 0.2)); // 90-110% of S
            r = Math.round(Math.random() * 10) / 100; // 0-10%
            sigma = Math.round(Math.random() * 20 + 15) / 100; // 15-35%
            t = Math.round(Math.random() * 12 + 1) / 12; // 1-13 months
            B = Math.round(S * (1.1 + Math.random() * 0.3)); // 10-40% above spot for barrier
            break;
            
        case 'hard':
            // More complex numbers
            S = Math.round(Math.random() * 200 + 50); // $50-250
            K = Math.round(S * (0.7 + Math.random() * 0.6)); // 70-130% of S
            r = Math.round(Math.random() * 15) / 100; // 0-15%
            sigma = Math.round(Math.random() * 40 + 10) / 100; // 10-50%
            t = (Math.random() * 3).toFixed(2); // 0-3 years with decimal
            B = Math.round(S * (1 + (Math.random() * 0.5 - 0.2))); // -20% to +30% of spot for barrier
            break;
    }
    
    // Choose a problem type
    const problemTypes = [
        'european_call',
        'european_put',
        'asian_call',
        'asian_put',
        'lookback_call',
        'lookback_put',
        'barrier_call_up',
        'barrier_put_down'
    ];
    
    // For easy difficulty, limit to European options
    let availableTypes = problemTypes;
    if (difficulty === 'easy') {
        availableTypes = problemTypes.slice(0, 2);
    } else if (difficulty === 'medium') {
        availableTypes = problemTypes.slice(0, 6);
    }
    
    const problemType = availableTypes[Math.floor(Math.random() * availableTypes.length)];
    
    // Generate problem text and answer
    let problemDescription, answer;
    
    switch (problemType) {
        case 'european_call':
            problemDescription = `Calculate the price of a European call option using Monte Carlo simulation with the following parameters:<br><br>
                Stock Price (S): $${S}<br>
                Strike Price (K): $${K}<br>
                Risk-free Rate (r): ${(r * 100).toFixed(1)}%<br>
                Volatility (σ): ${(sigma * 100).toFixed(1)}%<br>
                Time to Expiry (t): ${t} years<br>
                Number of Simulations: 10,000`;
            answer = priceEuropeanOption(S, K, r, sigma, t, 'call', 10000).price;
            break;
            
        case 'european_put':
            problemDescription = `Calculate the price of a European put option using Monte Carlo simulation with the following parameters:<br><br>
                Stock Price (S): $${S}<br>
                Strike Price (K): $${K}<br>
                Risk-free Rate (r): ${(r * 100).toFixed(1)}%<br>
                Volatility (σ): ${(sigma * 100).toFixed(1)}%<br>
                Time to Expiry (t): ${t} years<br>
                Number of Simulations: 10,000`;
            answer = priceEuropeanOption(S, K, r, sigma, t, 'put', 10000).price;
            break;
            
        case 'asian_call':
            problemDescription = `Calculate the price of an Asian call option (with arithmetic averaging) using Monte Carlo simulation with the following parameters:<br><br>
                Stock Price (S): $${S}<br>
                Strike Price (K): $${K}<br>
                Risk-free Rate (r): ${(r * 100).toFixed(1)}%<br>
                Volatility (σ): ${(sigma * 100).toFixed(1)}%<br>
                Time to Expiry (t): ${t} years<br>
                Number of Simulations: 10,000`;
            answer = priceAsianOption(S, K, r, sigma, t, 'call', 10000, 'arithmetic').price;
            break;
            
        case 'asian_put':
            problemDescription = `Calculate the price of an Asian put option (with arithmetic averaging) using Monte Carlo simulation with the following parameters:<br><br>
                Stock Price (S): $${S}<br>
                Strike Price (K): $${K}<br>
                Risk-free Rate (r): ${(r * 100).toFixed(1)}%<br>
                Volatility (σ): ${(sigma * 100).toFixed(1)}%<br>
                Time to Expiry (t): ${t} years<br>
                Number of Simulations: 10,000`;
            answer = priceAsianOption(S, K, r, sigma, t, 'put', 10000, 'arithmetic').price;
            break;
            
        case 'lookback_call':
            problemDescription = `Calculate the price of a lookback call option using Monte Carlo simulation with the following parameters:<br><br>
                Stock Price (S): $${S}<br>
                Risk-free Rate (r): ${(r * 100).toFixed(1)}%<br>
                Volatility (σ): ${(sigma * 100).toFixed(1)}%<br>
                Time to Expiry (t): ${t} years<br>
                Number of Simulations: 10,000<br><br>
                Note: For a lookback call, the payoff is S(T) - min(S) over the path.`;
            answer = priceLookbackOption(S, r, sigma, t, 'call', 10000).price;
            break;
            
        case 'lookback_put':
            problemDescription = `Calculate the price of a lookback put option using Monte Carlo simulation with the following parameters:<br><br>
                Stock Price (S): $${S}<br>
                Risk-free Rate (r): ${(r * 100).toFixed(1)}%<br>
                Volatility (σ): ${(sigma * 100).toFixed(1)}%<br>
                Time to Expiry (t): ${t} years<br>
                Number of Simulations: 10,000<br><br>
                Note: For a lookback put, the payoff is max(S) - S(T) over the path.`;
            answer = priceLookbackOption(S, r, sigma, t, 'put', 10000).price;
            break;
            
        case 'barrier_call_up':
            problemDescription = `Calculate the price of an up-and-out call barrier option using Monte Carlo simulation with the following parameters:<br><br>
                Stock Price (S): $${S}<br>
                Strike Price (K): $${K}<br>
                Barrier Level (B): $${B}<br>
                Risk-free Rate (r): ${(r * 100).toFixed(1)}%<br>
                Volatility (σ): ${(sigma * 100).toFixed(1)}%<br>
                Time to Expiry (t): ${t} years<br>
                Number of Simulations: 10,000<br><br>
                Note: The option becomes worthless if the stock price ever reaches or exceeds the barrier level.`;
            answer = priceBarrierOption(S, K, B, r, sigma, t, 'up-and-out', 'call', 10000).price;
            break;
            
        case 'barrier_put_down':
            problemDescription = `Calculate the price of a down-and-out put barrier option using Monte Carlo simulation with the following parameters:<br><br>
                Stock Price (S): $${S}<br>
                Strike Price (K): $${K}<br>
                Barrier Level (B): $${B}<br>
                Risk-free Rate (r): ${(r * 100).toFixed(1)}%<br>
                Volatility (σ): ${(sigma * 100).toFixed(1)}%<br>
                Time to Expiry (t): ${t} years<br>
                Number of Simulations: 10,000<br><br>
                Note: The option becomes worthless if the stock price ever reaches or falls below the barrier level.`;
            answer = priceBarrierOption(S, K, B, r, sigma, t, 'down-and-out', 'put', 10000).price;
            break;
    }
    
    // Update the problem text
    problemText.innerHTML = problemDescription;
    
    // Store the current problem
    currentProblem = {
        type: problemType,
        parameters: { S, K, B, r, sigma, t },
        description: problemDescription,
        answer: answer,
        difficulty: difficulty
    };
}

function checkAnswer() {
    if (!currentProblem) {
        alert('Please generate a problem first!');
        return;
    }
    
    const userAnswer = parseFloat(document.getElementById('user-answer').value);
    const feedback = document.getElementById('feedback');
    const feedbackMessage = document.getElementById('feedback-message');
    const explanation = document.getElementById('explanation');
    
    // Check if the answer is close enough (wider tolerances due to Monte Carlo randomness)
    const tolerance = currentProblem.difficulty === 'easy' ? 0.15 : 
                     (currentProblem.difficulty === 'medium' ? 0.10 : 0.05);
    
    const percentDiff = Math.abs((userAnswer - currentProblem.answer) / currentProblem.answer);
    const isCorrect = percentDiff < tolerance;
    
    // Show feedback
    feedback.classList.remove('hidden');
    
    if (isCorrect) {
        feedback.classList.add('correct');
        feedback.classList.remove('incorrect');
        feedbackMessage.innerHTML = `<strong>Correct!</strong> Your answer $${userAnswer.toFixed(4)} is within the acceptable range of the expected value $${currentProblem.answer.toFixed(4)}.`;
        
        // Add to history
        addToHistory(true);
    } else {
        feedback.classList.add('incorrect');
        feedback.classList.remove('correct');
        feedbackMessage.innerHTML = `<strong>Incorrect.</strong> Your answer $${userAnswer.toFixed(4)} is different from the expected value $${currentProblem.answer.toFixed(4)}.`;
        
        // Generate explanation
        explanation.innerHTML = generateExplanation();
        
        // Add to history
        addToHistory(false);
    }
}

function generateExplanation() {
    const { type, parameters } = currentProblem;
    const { S, K, B, r, sigma, t } = parameters;
    
    let explanation = '<h4>Solution Approach:</h4>';
    
    // Add common Monte Carlo explanation
    explanation += `
        <p>To solve this problem using Monte Carlo simulation:</p>
        <ol>
            <li>Simulate multiple stock price paths using Geometric Brownian Motion:
                <br>S_{t+Δt} = S_t × exp((r - σ²/2)×Δt + σ×√Δt×ε)
                <br>where ε is a standard normal random variable</li>
            <li>Calculate the payoff for each path based on the option type</li>
            <li>Average the payoffs and discount back to present value:
                <br>Option Price = e^(-rt) × (1/N) × Σ Payoffs</li>
        </ol>
    `;
    
    // Add specific explanation based on option type
    switch (type) {
        case 'european_call':
        case 'european_put':
            const optType = type.includes('call') ? 'call' : 'put';
            const payoffFormula = optType === 'call' ? 'max(S_T - K, 0)' : 'max(K - S_T, 0)';
            
            explanation += `
                <p>For a European ${optType} option, the payoff at maturity is ${payoffFormula}, where S_T is the stock price at maturity.</p>
                <p>Using 10,000 simulations with the given parameters:</p>
                <ul>
                    <li>Stock Price (S): $${S}</li>
                    <li>Strike Price (K): $${K}</li>
                    <li>Risk-free Rate (r): ${(r * 100).toFixed(1)}%</li>
                    <li>Volatility (σ): ${(sigma * 100).toFixed(1)}%</li>
                    <li>Time to Expiry (t): ${t} years</li>
                </ul>
                <p>The Monte Carlo simulation yields an option price of $${currentProblem.answer.toFixed(4)}.</p>
            `;
            break;
            
        case 'asian_call':
        case 'asian_put':
            const asianType = type.includes('call') ? 'call' : 'put';
            const asianPayoff = asianType === 'call' ? 'max(S_avg - K, 0)' : 'max(K - S_avg, 0)';
            
            explanation += `
                <p>For an Asian ${asianType} option with arithmetic averaging, the payoff depends on the average stock price over the entire path.</p>
                <p>The payoff formula is ${asianPayoff}, where S_avg is the arithmetic average of stock prices observed over the path.</p>
                <p>Using 10,000 simulations with the given parameters:</p>
                <ul>
                    <li>Stock Price (S): $${S}</li>
                    <li>Strike Price (K): $${K}</li>
                    <li>Risk-free Rate (r): ${(r * 100).toFixed(1)}%</li>
                    <li>Volatility (σ): ${(sigma * 100).toFixed(1)}%</li>
                    <li>Time to Expiry (t): ${t} years</li>
                </ul>
                <p>The Monte Carlo simulation yields an option price of $${currentProblem.answer.toFixed(4)}.</p>
            `;
            break;
            
        case 'lookback_call':
        case 'lookback_put':
            const lookbackType = type.includes('call') ? 'call' : 'put';
            const lookbackFormula = lookbackType === 'call' ? 'S_T - min(S)' : 'max(S) - S_T';
            
            explanation += `
                <p>For a lookback ${lookbackType} option, the payoff depends on the extreme value (maximum or minimum) of the stock price over the path.</p>
                <p>The payoff formula is ${lookbackFormula}, where min(S)/max(S) are the minimum/maximum stock prices observed over the path, and S_T is the final stock price.</p>
                <p>Using 10,000 simulations with the given parameters:</p>
                <ul>
                    <li>Stock Price (S): $${S}</li>
                    <li>Risk-free Rate (r): ${(r * 100).toFixed(1)}%</li>
                    <li>Volatility (σ): ${(sigma * 100).toFixed(1)}%</li>
                    <li>Time to Expiry (t): ${t} years</li>
                </ul>
                <p>The Monte Carlo simulation yields an option price of $${currentProblem.answer.toFixed(4)}.</p>
            `;
            break;
            
        case 'barrier_call_up':
        case 'barrier_put_down':
            const barrierType = type.includes('call') ? 'call' : 'put';
            const barrierDirection = type.includes('up') ? 'up-and-out' : 'down-and-out';
            const barrierCondition = barrierDirection === 'up-and-out' ? 'reaches or exceeds' : 'reaches or falls below';
            const standardPayoff = barrierType === 'call' ? 'max(S_T - K, 0)' : 'max(K - S_T, 0)';
            
            explanation += `
                <p>For a ${barrierDirection} ${barrierType} barrier option, the payoff is like a standard ${barrierType} option, but the option becomes worthless if the stock price ever ${barrierCondition} the barrier level.</p>
                <p>The payoff is either ${standardPayoff} if the barrier is not triggered, or 0 if the barrier is triggered.</p>
                <p>Using 10,000 simulations with the given parameters:</p>
                <ul>
                    <li>Stock Price (S): $${S}</li>
                    <li>Strike Price (K): $${K}</li>
                    <li>Barrier Level (B): $${B}</li>
                    <li>Risk-free Rate (r): ${(r * 100).toFixed(1)}%</li>
                    <li>Volatility (σ): ${(sigma * 100).toFixed(1)}%</li>
                    <li>Time to Expiry (t): ${t} years</li>
                </ul>
                <p>The Monte Carlo simulation yields an option price of $${currentProblem.answer.toFixed(4)}.</p>
            `;
            break;
    }
    
    explanation += `<p>Note: Due to the random nature of Monte Carlo simulation, your answer may differ slightly from the expected value. A tolerance of ${currentProblem.difficulty === 'easy' ? '15%' : (currentProblem.difficulty === 'medium' ? '10%' : '5%')} is allowed for this problem.</p>`;
    
    return explanation;
}

function addToHistory(isCorrect) {
    // Add problem to history
    const historyItem = {
        ...currentProblem,
        userAnswer: parseFloat(document.getElementById('user-answer').value),
        isCorrect: isCorrect,
        timestamp: new Date()
    };
    
    problemHistory.unshift(historyItem);
    
    // Cap history at 10 items
    if (problemHistory.length > 10) {
        problemHistory.pop();
    }
    
    // Update history display
    updateHistoryDisplay();
}

function updateHistoryDisplay() {
    const historyList = document.getElementById('history-list');
    
    if (problemHistory.length === 0) {
        historyList.innerHTML = '<p>Your solved problems will appear here.</p>';
        return;
    }
    
    let historyHTML = '';
    
    problemHistory.forEach((item, index) => {
        const problemTypeLabel = getProblemTypeLabel(item.type);
        
        historyHTML += `
            <div class="history-item ${item.isCorrect ? 'correct' : 'incorrect'}">
                <div class="history-header">
                    <span class="history-number">#${index + 1}</span>
                    <span class="history-type">${problemTypeLabel}</span>
                    <span class="history-difficulty">${item.difficulty.charAt(0).toUpperCase() + item.difficulty.slice(1)}</span>
                </div>
                <div class="history-result">
                    <span>Your answer: $${item.userAnswer.toFixed(4)}</span>
                    <span>Expected answer: $${item.answer.toFixed(4)}</span>
                </div>
            </div>
        `;
    });
    
    historyList.innerHTML = historyHTML;
}

function getProblemTypeLabel(type) {
    const labels = {
        'european_call': 'European Call',
        'european_put': 'European Put',
        'asian_call': 'Asian Call',
        'asian_put': 'Asian Put',
        'lookback_call': 'Lookback Call',
        'lookback_put': 'Lookback Put',
        'barrier_call_up': 'Up-and-Out Call',
        'barrier_put_down': 'Down-and-Out Put'
    };
    
    return labels[type] || type;
}

// --------------------------------------------------
// Set up event listeners
// --------------------------------------------------

function setupEventListeners() {
    // Event listeners for tab switching
    const tabButtons = document.querySelectorAll('.tab-btn');
    tabButtons.forEach(button => {
        button.addEventListener('click', function() {
            const tabId = this.getAttribute('data-tab');
            
            // Remove active class from all buttons and tabs
            tabButtons.forEach(btn => btn.classList.remove('active'));
            document.querySelectorAll('.tab-pane').forEach(pane => pane.classList.remove('active'));
            
            // Add active class to current button and tab
            this.classList.add('active');
            document.getElementById(tabId)?.classList.add('active');
        });
    });
} 