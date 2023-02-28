#pragma once

#define CUPCT1_M 5
#define CUPCT1_P 3
#define CUPCT1_W 2
#define CUPCT1_CORRSIZE 40
#define CUPCT1_ADJSIZE 59

/**
 * This is a row-major matrix with the following structure:
 * m + p rows, w + p cols,
 * {
 * 0.47, 0.39, 0.74, 0.53, 0.62,
 * 0.64, 0.51, 0.35, 0.4,  0.15,
 * 0.80, 0.15, 0.98, 0.21, 0.94,
 * 0.50, 0.00, 0.33, 0.49, 0.25, // here only one forward entry
 * 0.00, 0.00, 0.28, 0.38, 0.24, // here no forward entries
 * 0.00, 0.00, 0.00, 0.21, 0.06, // the last three lines are the phenotype corrs in the lower right corner
 * 0.00, 0.00, 0.21, 0.00, 0.99,
 * 0.00, 0.00, 0.06, 0.99, 0.00
 * }
 */
const float cupct1_c[40] = {0.47, 0.39, 0.74, 0.53, 0.62, 0.64, 0.51, 0.35, 0.4, 0.15, 0.8,
                            0.15, 0.98, 0.21, 0.94, 0.5, 0., 0.33, 0.49, 0.25, 0., 0.,
                            0.28, 0.38, 0.24, 0., 0., 0., 0.21, 0.06, 0., 0., 0.21,
                            0., 0.99, 0., 0., 0.06, 0.99, 0.};

/**
 * This is a row-major compound matrix that consists of two parts.
 * First part:
 * m rows, 2 * w + p cols
 * {
 * 0, 0, 0, 0, 1, 0, 0,
 * 0, 0, 1, 0, 0, 0, 0,
 * 0, 1, 1, 0, 1, 0, 1,
 * 0, 1, 0, 0, 0, 0, 0,
 * 0, 0, 0, 0, 0, 0, 0,
 * }
 * Second part:
 * p rows, m + p cols
 * {
 * 1, 0, 1, 0, 0, 0, 0, 0,
 * 0, 0, 0, 0, 0, 0, 0, 1,
 * 0, 0, 1, 0, 0, 0, 1, 0
 * }
 */
const int cupct1_g[59] = {0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 1, 0,
                          1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0,
                          0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0};

#define CUPCT2_M 200
#define CUPCT2_P 3
#define CUPCT2_W 10
#define CUPCT2_CORRSIZE 2639
#define CUPCT2_ADJSIZE 5209

const float cupct2_c[CUPCT2_CORRSIZE] = {0.42, 0.87, 0.95, 0.28, 0.69, 0.83, 0.2, 0.17, 0.13, 0.25, 0.19, 0.84, 0.95, 0.32, 0.95, 0.08, 0.97, 0.56, 0.14, 0.13, 0.68, 0.26, 0.69, 0.32, 0.14, 0.1, 0.01, 0.65, 0.14, 0.2, 0.4, 0.7, 0.41, 0.17, 0.51, 0.62, 0.22, 0.37, 0.26, 0.73, 0.24, 0.56, 0.13, 0.03, 0.02, 0.82, 0.17, 0.35, 0.83, 0.38, 0.81, 0.21, 0.05, 0.75, 0.19, 0.4, 0.56, 0.5, 0.67, 0.06, 0.31, 0.1, 0.32, 0.48, 0.24, 0.62, 0.22, 0.65, 0.2, 0.13, 0.49, 0.68, 0.93, 0.29, 0.78, 0.65, 0.56, 0.46, 0.98, 0.94, 0.77, 0.56, 0.27, 0.49, 0.22, 0.54, 0.15, 0.91, 0.65, 0.95, 0.84, 0.79, 0.22, 0.06, 0.4, 0.28, 0.38, 0.17, 0.74, 0.18, 0.28, 0.84, 0.23, 0.74, 0.58, 0.69, 0.1, 0.71, 0.24, 0.42, 0.79, 0.79, 0.35, 0.83, 0.55, 0.64, 0.68, 0.5, 0.51, 0.44, 0.54, 0.46, 0.69, 0.05, 0.85, 0.75, 0.64, 0.86, 0.77, 0.56, 0.03, 0.43, 0.64, 0.53, 0.76, 0.18, 0.85, 0.61, 0.87, 0.15, 0.95, 0.01, 0.52, 0.83, 0.22, 0.01, 0.19, 0.1, 0.68, 0.01, 0.46, 0.47, 0.83, 0.96, 0.61, 0.43, 0.9, 0.2, 0.43, 0.79, 0.17, 0.85, 0.77, 0.24, 0.47, 0.98, 0.72, 0.3, 0.89, 0.56, 0.39, 0.39, 0.25, 0.17, 0.68, 0.37, 0.17, 0.95, 0.77, 0.16, 0.61, 0.78, 0.8, 0.46, 0.99, 0.87, 0.47, 0.89, 0.52, 0.59, 0.95, 0.22, 0.68, 0.2, 0.48, 0.29, 1.0, 0.23, 0.72, 0.94, 0.07, 0.85, 0.26, 0.49, 0.63, 0.66, 0.96, 0.11, 0.53, 0.64, 0.65, 0.54, 0.35, 0.09, 0.23, 0.41, 0.41, 0.78, 0.66, 0.1, 0.84, 0.97, 0.98, 0.35, 0.47, 0.96, 0.96, 0.9, 0.79, 0.27, 0.98, 0.43, 0.43, 0.07, 0.35, 0.86, 0.65, 0.88, 0.92, 0.38, 0.55, 0.98, 0.83, 0.61, 0.71, 0.13, 0.88, 0.4, 0.04, 0.55, 0.23, 0.2, 0.93, 0.72, 0.16, 0.05, 0.9, 0.11, 0.7, 0.6, 0.05, 0.62, 0.44, 0.14, 0.96, 0.14, 0.81, 0.4, 0.98, 0.67, 0.88, 0.44, 0.9, 0.12, 0.73, 0.64, 0.75, 0.81, 0.31, 0.18, 0.62, 0.32, 0.05, 0.28, 0.55, 0.4, 0.06, 0.29, 0.55, 0.34, 0.71, 0.84, 0.17, 0.77, 0.95, 0.69, 0.23, 0.86, 0.27, 0.08, 0.27, 0.4, 0.35, 0.6, 0.35, 0.45, 0.43, 0.07, 0.92, 0.74, 0.33, 0.1, 0.44, 0.57, 0.34, 0.78, 0.31, 0.73, 0.8, 0.31, 0.06, 0.47, 0.42, 0.82, 0.05, 0.8, 0.54, 0.04, 0.46, 0.22, 0.72, 0.28, 0.65, 0.54, 0.9, 0.16, 0.71, 0.5, 0.62, 0.7, 0.83, 0.47, 0.22, 0.9, 0.6, 0.61, 0.88, 0.97, 0.99, 0.58, 0.4, 0.08, 0.38, 0.11, 0.38, 0.12, 0.03, 0.35, 0.08, 0.67, 0.64, 0.05, 0.22, 0.95, 0.69, 0.33, 0.18, 0.57, 0.56, 0.55, 0.23, 0.95, 0.87, 0.47, 0.53, 0.52, 0.16, 0.42, 0.97, 0.09, 0.56, 0.91, 0.34, 0.74, 0.8, 0.89, 0.97, 0.47, 0.59, 0.47, 0.07, 0.71, 0.46, 0.17, 0.87, 0.54, 0.22, 0.9, 0.67, 0.39, 0.59, 0.91, 0.72, 0.02, 0.22, 0.29, 0.74, 0.19, 0.48, 0.84, 0.04, 0.28, 0.84, 0.85, 0.84, 0.52, 0.76, 0.59, 0.11, 0.89, 0.94, 0.53, 0.4, 0.27, 0.72, 0.37, 0.08, 0.45, 0.17, 0.1, 0.87, 0.28, 0.9, 0.68, 0.12, 0.09, 0.83, 0.81, 0.65, 0.19, 0.85, 0.57, 0.33, 0.33, 0.35, 0.59, 0.27, 0.59, 0.32, 0.98, 0.16, 0.03, 0.06, 0.63, 0.94, 0.97, 0.92, 0.09, 0.97, 0.68, 0.5, 0.37, 0.16, 0.55, 0.68, 0.31, 0.1, 0.64, 0.26, 0.39, 0.43, 0.28, 0.93, 0.49, 0.24, 0.75, 0.22, 0.61, 0.21, 0.84, 0.98, 0.31, 0.86, 0.84, 0.35, 0.65, 0.76, 0.83, 0.23, 0.61, 0.78, 0.84, 0.99, 0.57, 0.9, 0.23, 0.44, 0.59, 0.18, 0.83, 0.16, 0.1, 0.21, 0.23, 0.7, 0.63, 0.66, 0.57, 0.22, 0.02, 0.36, 0.81, 0.93, 0.27, 0.8, 0.45, 0.71, 0.36, 0.87, 0.76, 0.5, 0.48, 0.24, 0.43, 0.42, 0.37, 0.91, 0.0, 0.2, 0.9, 0.11, 0.82, 0.96, 0.2, 0.48, 0.5, 0.37, 0.35, 0.41, 0.21, 0.88, 0.91, 0.3, 0.43, 0.56, 0.29, 0.17, 0.19, 0.86, 0.1, 0.19, 0.99, 0.35, 0.47, 0.23, 0.32, 0.67, 0.99, 0.02, 0.6, 0.51, 0.76, 0.5, 0.28, 0.56, 0.59, 0.61, 0.85, 0.83, 0.8, 0.48, 0.12, 0.79, 0.14, 0.04, 0.26, 0.03, 0.96, 0.24, 0.59, 0.43, 0.66, 0.77, 0.1, 0.24, 0.81, 0.86, 0.54, 0.49, 0.81, 0.05, 0.93, 0.25, 0.88, 0.42, 0.86, 0.3, 0.88, 0.13, 0.79, 0.43, 0.06, 0.62, 0.89, 0.68, 0.06, 0.31, 0.49, 0.1, 0.39, 0.91, 0.84, 0.08, 0.71, 0.97, 0.68, 0.49, 0.1, 0.43, 0.63, 0.53, 0.95, 0.3, 0.72, 0.71, 0.87, 0.5, 0.53, 0.04, 0.99, 0.52, 0.34, 0.3, 0.07, 0.93, 0.09, 0.62, 0.68, 0.38, 0.42, 0.06, 0.66, 0.44, 0.34, 0.21, 0.56, 0.92, 0.75, 0.56, 0.43, 0.24, 0.13, 0.68, 0.3, 0.32, 0.63, 0.05, 0.78, 0.98, 0.16, 0.9, 0.44, 0.97, 0.98, 0.72, 0.49, 0.98, 0.59, 0.52, 0.36, 0.24, 0.39, 0.08, 0.67, 0.69, 0.85, 0.93, 0.56, 0.74, 0.54, 0.05, 0.17, 0.04, 0.57, 0.91, 0.01, 0.86, 0.24, 0.24, 0.96, 0.54, 0.72, 0.05, 0.33, 0.1, 0.98, 0.19, 0.26, 0.73, 0.24, 0.08, 0.56, 0.09, 0.79, 0.27, 0.95, 0.75, 0.05, 0.71, 0.21, 0.81, 0.62, 0.92, 0.23, 0.86, 0.19, 0.13, 0.97, 0.75, 0.22, 0.61, 0.66, 0.47, 0.31, 0.06, 0.13, 0.48, 0.62, 0.84, 0.43, 0.85, 0.76, 0.13, 0.04, 0.83, 0.08, 0.86, 0.24, 0.41, 0.81, 0.19, 0.96, 0.7, 0.23, 0.78, 0.47, 0.35, 0.35, 0.1, 0.44, 0.54, 0.45, 0.37, 0.65, 0.71, 0.66, 0.42, 0.97, 0.84, 0.26, 0.28, 0.67, 0.82, 0.59, 0.51, 0.6, 0.77, 0.81, 0.73, 0.86, 0.03, 0.9, 0.56, 0.54, 0.58, 0.09, 0.03, 0.77, 0.67, 0.61, 0.69, 0.48, 0.58, 0.33, 0.52, 0.54, 0.96, 0.13, 0.83, 0.95, 0.83, 0.81, 0.51, 0.39, 0.41, 0.3, 0.52, 0.26, 0.36, 0.52, 0.02, 0.78, 0.54, 0.58, 0.43, 0.58, 0.49, 0.36, 0.96, 0.96, 0.28, 0.91, 0.84, 0.59, 0.64, 0.83, 0.76, 0.07, 0.7, 0.97, 0.59, 0.73, 0.05, 0.75, 0.79, 0.75, 0.61, 0.23, 0.47, 0.59, 0.81, 0.81, 0.87, 0.46, 0.78, 0.07, 0.07, 0.53, 0.58, 0.53, 0.47, 0.38, 0.95, 0.84, 0.11, 0.91, 0.63, 0.7, 0.43, 0.3, 0.89, 0.23, 0.72, 0.78, 0.9, 0.73, 0.16, 0.23, 0.72, 0.89, 0.45, 0.72, 0.74, 0.98, 0.73, 0.2, 0.31, 0.58, 0.95, 0.98, 0.1, 0.28, 0.17, 0.93, 0.34, 0.95, 0.77, 0.27, 0.42, 0.32, 0.58, 0.93, 0.45, 0.44, 0.55, 0.22, 0.64, 0.89, 0.46, 0.84, 0.12, 0.25, 0.38, 0.44, 0.99, 0.89, 0.75, 0.23, 0.36, 0.54, 0.02, 0.52, 0.4, 0.02, 0.3, 0.2, 0.9, 0.35, 0.08, 0.14, 0.97, 0.93, 0.21, 0.14, 0.03, 0.78, 0.08, 0.07, 0.46, 0.24, 0.23, 0.3, 0.69, 0.44, 0.75, 0.35, 0.55, 0.72, 0.43, 0.54, 0.91, 0.96, 0.84, 0.59, 0.81, 0.65, 0.33, 0.57, 0.78, 0.56, 0.48, 0.12, 0.8, 0.14, 0.95, 0.65, 0.93, 0.88, 0.65, 0.94, 0.31, 0.17, 0.94, 0.36, 0.2, 0.49, 0.95, 0.13, 0.61, 0.71, 0.69, 0.85, 0.04, 0.3, 0.33, 0.33, 0.12, 0.63, 0.1, 0.09, 0.85, 0.11, 0.66, 0.42, 0.74, 0.8, 0.83, 0.55, 0.16, 0.01, 0.0, 0.77, 0.21, 0.87, 0.95, 0.92, 0.59, 0.1, 0.37, 0.82, 0.37, 0.81, 0.46, 0.4, 0.27, 0.06, 0.21, 0.62, 0.91, 0.02, 0.88, 0.36, 0.96, 0.82, 0.71, 0.61, 0.62, 0.76, 0.19, 0.63, 0.36, 0.43, 0.52, 0.62, 0.16, 0.32, 0.4, 0.29, 0.52, 0.77, 0.87, 0.45, 0.43, 0.07, 0.68, 0.21, 0.36, 0.61, 0.02, 0.85, 0.72, 0.59, 0.53, 0.7, 0.47, 0.04, 0.88, 0.97, 0.49, 0.22, 0.46, 0.79, 0.32, 0.48, 0.5, 0.53, 0.23, 0.98, 0.18, 0.23, 0.8, 0.62, 0.59, 0.23, 0.49, 0.79, 0.49, 0.05, 0.06, 0.99, 0.54, 0.88, 0.73, 0.29, 0.37, 0.43, 0.36, 0.1, 0.1, 0.96, 0.35, 0.03, 0.29, 0.87, 0.23, 0.87, 0.46, 0.38, 0.9, 0.79, 0.14, 0.03, 0.06, 0.38, 0.86, 0.07, 0.91, 0.87, 1.0, 0.55, 0.66, 0.56, 0.74, 0.05, 0.21, 0.91, 0.07, 0.75, 0.01, 0.79, 0.95, 0.05, 0.37, 0.54, 0.06, 0.04, 0.76, 0.27, 0.09, 0.31, 0.08, 0.26, 0.06, 0.35, 0.13, 0.01, 0.36, 0.52, 0.06, 0.83, 0.53, 0.74, 0.38, 0.37, 0.63, 0.42, 0.04, 0.24, 0.86, 0.33, 0.53, 0.52, 0.02, 0.12, 0.75, 0.37, 0.68, 0.39, 0.63, 0.24, 0.84, 1.0, 0.86, 0.03, 0.03, 0.46, 0.98, 0.69, 0.09, 0.18, 0.17, 0.36, 0.79, 0.72, 0.62, 0.75, 0.16, 0.5, 0.56, 0.4, 0.94, 0.84, 0.15, 0.58, 0.98, 0.7, 0.69, 0.74, 0.59, 0.07, 0.88, 0.44, 0.25, 0.11, 0.83, 0.64, 0.26, 0.91, 0.56, 0.19, 0.3, 0.88, 0.94, 0.98, 0.52, 0.44, 0.88, 0.99, 0.51, 0.19, 0.83, 0.91, 0.72, 0.92, 0.23, 0.72, 0.89, 0.3, 0.08, 0.06, 0.54, 0.48, 0.66, 0.69, 0.12, 0.29, 0.04, 0.06, 0.57, 0.98, 0.36, 0.41, 0.9, 0.04, 0.95, 0.69, 0.21, 0.37, 0.76, 0.16, 0.43, 0.68, 0.8, 0.17, 0.14, 0.67, 0.92, 0.97, 0.05, 0.69, 0.13, 0.91, 0.03, 0.73, 0.8, 0.54, 0.19, 0.46, 0.16, 0.28, 0.37, 0.03, 0.49, 0.49, 0.32, 0.98, 0.86, 0.42, 0.44, 0.48, 0.01, 0.02, 0.96, 0.95, 0.24, 0.26, 1.0, 0.93, 0.26, 0.14, 0.18, 0.62, 0.67, 0.71, 0.31, 0.96, 0.43, 0.21, 0.92, 0.71, 0.01, 0.91, 0.73, 0.26, 0.12, 0.63, 0.82, 0.34, 0.41, 0.82, 0.86, 0.22, 0.06, 0.54, 0.93, 0.6, 0.96, 0.0, 0.85, 0.7, 0.23, 0.9, 0.22, 0.77, 0.51, 0.29, 0.8, 0.62, 0.88, 0.47, 0.9, 0.52, 0.66, 0.18, 0.38, 0.98, 0.17, 0.64, 0.94, 0.69, 0.44, 0.28, 0.02, 0.47, 0.49, 0.39, 0.81, 1.0, 0.16, 0.07, 0.84, 0.63, 0.19, 0.86, 0.98, 0.7, 0.28, 0.99, 0.01, 0.58, 0.59, 0.79, 0.79, 0.71, 0.37, 0.1, 0.24, 0.56, 0.27, 0.38, 0.71, 0.36, 0.53, 0.37, 0.96, 0.4, 0.87, 0.83, 0.45, 0.74, 0.24, 0.13, 0.26, 0.33, 0.88, 0.65, 0.74, 0.82, 0.52, 0.69, 0.38, 0.7, 0.03, 0.89, 0.73, 0.32, 0.87, 0.54, 0.8, 0.92, 0.02, 0.9, 0.47, 0.75, 0.15, 0.03, 0.54, 0.68, 0.48, 0.25, 0.95, 0.63, 0.13, 0.48, 0.5, 0.24, 0.78, 0.54, 0.25, 0.36, 0.66, 0.14, 0.16, 0.4, 0.24, 0.39, 0.88, 0.81, 0.57, 0.72, 0.84, 0.12, 0.06, 0.42, 0.16, 0.16, 0.77, 0.57, 0.5, 0.16, 0.98, 0.18, 0.46, 0.28, 0.38, 0.29, 0.04, 0.4, 0.43, 0.18, 0.88, 0.97, 0.7, 0.72, 0.18, 0.57, 0.81, 0.29, 0.1, 0.34, 0.08, 0.77, 0.44, 0.66, 0.84, 0.35, 0.3, 0.39, 0.43, 0.07, 0.62, 0.93, 0.66, 0.91, 0.45, 0.84, 0.09, 0.25, 0.13, 0.5, 0.81, 0.98, 0.75, 0.79, 0.51, 0.46, 0.44, 0.63, 0.33, 0.07, 0.21, 0.57, 0.65, 0.56, 0.43, 0.47, 0.62, 0.77, 0.24, 0.53, 0.37, 0.55, 0.9, 0.32, 0.4, 0.04, 0.43, 0.42, 0.19, 0.31, 0.31, 0.19, 0.9, 0.93, 0.83, 0.53, 0.05, 0.37, 0.47, 0.55, 0.73, 0.2, 0.27, 0.51, 0.76, 0.03, 0.74, 0.8, 0.52, 0.12, 0.93, 0.35, 0.67, 0.78, 0.91, 0.49, 0.51, 0.43, 0.05, 0.73, 0.19, 0.45, 0.68, 0.13, 0.42, 0.01, 0.4, 0.5, 0.42, 0.67, 0.67, 0.77, 0.39, 0.92, 0.14, 0.14, 0.05, 0.69, 0.21, 0.55, 0.24, 0.44, 0.47, 0.38, 0.61, 0.61, 0.06, 0.31, 0.64, 0.51, 0.04, 0.85, 0.09, 0.84, 0.96, 0.48, 0.56, 0.34, 0.74, 0.77, 0.42, 0.11, 0.52, 0.22, 0.86, 0.82, 0.76, 0.64, 0.13, 0.78, 0.81, 0.49, 0.34, 0.25, 0.77, 0.75, 0.23, 0.13, 0.5, 0.57, 0.98, 0.14, 0.51, 0.93, 0.27, 0.95, 0.37, 0.17, 0.73, 0.81, 0.42, 0.14, 0.83, 0.41, 0.96, 0.15, 0.47, 0.22, 0.36, 0.97, 0.65, 0.15, 0.39, 0.29, 0.97, 0.12, 0.89, 0.48, 0.17, 0.23, 0.12, 0.62, 0.84, 0.87, 0.11, 0.91, 0.66, 0.7, 0.58, 0.95, 0.85, 0.23, 0.1, 0.42, 0.76, 0.12, 0.37, 0.71, 0.55, 0.7, 0.94, 0.53, 0.83, 0.38, 0.73, 0.21, 0.19, 0.66, 0.08, 0.75, 0.43, 0.6, 0.63, 0.77, 0.21, 0.31, 0.73, 0.57, 0.81, 0.39, 0.09, 0.44, 0.48, 0.35, 0.97, 0.13, 0.84, 0.06, 0.65, 0.45, 0.63, 0.43, 0.43, 0.49, 0.01, 0.41, 0.75, 0.87, 0.53, 0.43, 0.24, 0.69, 0.5, 0.76, 0.78, 0.6, 0.78, 0.43, 0.98, 0.29, 0.91, 0.3, 0.88, 0.87, 0.85, 0.54, 0.23, 0.39, 0.14, 0.47, 0.82, 0.72, 0.64, 0.26, 0.62, 0.11, 0.66, 0.29, 0.62, 0.19, 0.82, 0.1, 0.02, 0.24, 0.32, 0.09, 0.8, 0.39, 0.53, 0.43, 0.54, 0.45, 0.84, 0.44, 0.38, 0.29, 0.26, 0.95, 0.68, 0.17, 0.87, 0.7, 0.49, 0.0, 0.76, 0.02, 0.93, 0.88, 0.04, 0.87, 0.06, 0.38, 0.62, 0.98, 0.12, 0.63, 0.26, 0.33, 0.21, 0.45, 0.59, 0.49, 0.76, 0.54, 0.2, 0.45, 0.52, 0.03, 0.5, 0.12, 0.64, 0.53, 0.03, 0.75, 0.83, 0.16, 0.19, 0.12, 0.08, 0.5, 0.51, 0.53, 0.12, 0.56, 0.09, 0.62, 0.14, 0.47, 0.47, 0.92, 0.75, 0.24, 0.29, 0.66, 0.88, 0.02, 0.28, 0.81, 0.97, 0.91, 0.44, 0.28, 0.4, 0.29, 0.93, 0.15, 0.06, 0.18, 0.75, 0.21, 0.26, 0.8, 0.85, 0.66, 0.9, 0.02, 0.01, 0.65, 0.98, 0.26, 0.38, 0.6, 0.29, 0.04, 0.75, 0.45, 0.4, 0.43, 0.49, 0.26, 0.63, 0.02, 0.67, 0.57, 1.0, 0.56, 0.48, 0.41, 0.14, 0.52, 0.2, 0.56, 0.08, 0.35, 0.99, 0.79, 0.89, 0.69, 0.48, 0.98, 0.6, 0.14, 0.26, 0.44, 0.42, 0.39, 0.32, 0.02, 0.72, 0.66, 0.05, 0.89, 0.08, 0.13, 0.5, 0.91, 0.7, 0.43, 0.39, 0.49, 0.01, 0.34, 0.59, 0.13, 0.07, 0.09, 0.08, 0.87, 0.57, 0.22, 0.64, 0.4, 0.97, 0.43, 0.69, 0.07, 0.41, 0.82, 0.97, 0.78, 0.76, 0.11, 0.49, 0.76, 0.62, 0.77, 0.5, 0.02, 0.4, 0.67, 0.84, 0.13, 0.58, 0.55, 0.34, 0.69, 0.1, 0.85, 0.2, 0.24, 0.06, 0.31, 0.82, 0.62, 0.5, 0.83, 0.3, 0.79, 0.69, 0.83, 0.13, 0.82, 0.56, 0.43, 0.53, 0.3, 0.21, 0.9, 0.52, 0.88, 0.23, 0.78, 0.06, 0.45, 0.95, 0.18, 0.14, 0.95, 0.19, 0.36, 0.57, 0.53, 0.89, 0.16, 0.86, 0.4, 1.0, 0.65, 0.55, 0.82, 0.5, 0.37, 0.31, 0.22, 0.49, 0.16, 0.4, 0.94, 0.5, 0.98, 0.75, 0.57, 0.33, 0.1, 0.59, 0.56, 0.25, 0.89, 0.48, 0.28, 0.42, 0.69, 0.08, 0.85, 0.48, 0.04, 0.97, 0.83, 0.8, 0.81, 0.02, 0.57, 0.61, 0.55, 0.17, 0.1, 0.99, 0.4, 0.15, 0.47, 0.49, 0.68, 0.04, 0.12, 0.29, 0.78, 0.13, 0.28, 0.51, 0.57, 0.19, 0.58, 0.19, 0.38, 0.91, 0.17, 0.94, 0.7, 0.45, 0.32, 0.99, 0.47, 0.54, 0.9, 0.2, 0.63, 0.57, 0.06, 0.58, 0.77, 0.77, 0.8, 0.78, 0.48, 0.31, 0.92, 0.89, 0.15, 0.78, 0.08, 0.03, 0.79, 0.75, 0.53, 0.76, 0.38, 0.89, 0.73, 0.47, 0.79, 0.57, 0.58, 0.99, 0.5, 0.83, 0.22, 0.14, 0.28, 0.85, 0.52, 0.78, 0.25, 0.12, 0.83, 0.02, 0.61, 0.11, 0.65, 0.91, 0.77, 0.9, 0.98, 0.32, 0.13, 0.5, 0.74, 0.13, 0.56, 0.6, 0.32, 0.88, 0.83, 0.77, 0.46, 0.9, 0.26, 0.48, 0.27, 0.84, 0.65, 0.65, 0.76, 0.17, 0.66, 0.83, 0.2, 0.64, 0.93, 0.59, 0.41, 0.61, 0.96, 0.52, 0.19, 0.02, 0.43, 0.58, 0.56, 0.93, 0.94, 0.69, 0.46, 0.81, 0.56, 0.34, 0.17, 0.29, 0.71, 0.69, 0.66, 0.1, 0.01, 0.4, 0.48, 0.54, 0.22, 0.75, 0.53, 0.02, 0.93, 0.98, 0.55, 0.93, 0.8, 0.31, 0.78, 0.25, 0.33, 0.51, 0.17, 0.17, 0.29, 0.43, 0.43, 0.99, 0.12, 0.59, 0.03, 0.5, 0.52, 0.57, 0.76, 0.77, 0.61, 0.91, 0.29, 0.11, 0.47, 0.22, 0.12, 0.06, 0.29, 0.38, 0.69, 0.57, 0.01, 0.99, 0.86, 0.28, 0.69, 0.7, 0.6, 0.52, 0.39, 0.89, 0.08, 0.35, 0.45, 0.96, 0.41, 0.15, 0.08, 0.5, 0.35, 0.68, 0.27, 0.38, 0.57, 0.68, 0.52, 0.11, 0.54, 0.9, 0.15, 0.69, 0.68, 0.49, 0.31, 0.14, 0.53, 0.01, 0.03, 0.11, 0.85, 0.66, 0.54, 0.0, 0.74, 0.13, 0.31, 0.53, 0.83, 0.59, 0.89, 0.06, 0.88, 0.31, 0.34, 0.04, 0.48, 0.27, 0.99, 0.38, 0.56, 0.11, 0.59, 0.7, 0.54, 0.36, 0.41, 0.02, 0.34, 0.7, 0.6, 0.49, 0.91, 0.62, 0.18, 0.57, 0.4, 0.55, 0.77, 0.94, 0.36, 0.87, 0.14, 0.12, 0.74, 0.88, 0.6, 0.69, 0.45, 0.84, 0.9, 0.41, 0.8, 0.87, 0.95, 0.89, 0.23, 0.45, 0.87, 0.45, 0.79, 0.79, 0.68, 0.18, 0.85, 0.33, 0.35, 0.45, 0.48, 0.25, 0.39, 0.69, 0.28, 0.6, 0.43, 0.33, 0.17, 0.34, 0.1, 0.06, 0.55, 0.58, 0.67, 0.31, 0.64, 0.82, 0.82, 0.4, 0.39, 0.27, 0.23, 0.26, 0.9, 0.55, 0.53, 0.31, 0.64, 0.37, 0.11, 0.27, 0.37, 0.73, 0.68, 0.8, 0.99, 0.38, 0.29, 0.43, 0.33, 0.8, 0.44, 0.56, 0.15, 0.84, 0.4, 0.06, 0.5, 0.01, 0.72, 0.14, 0.33, 0.4, 0.54, 0.63, 0.69, 0.56, 0.25, 0.4, 0.74, 0.45, 0.64, 0.33, 0.89, 0.27, 0.84, 0.37, 0.63, 0.8, 0.61, 0.58, 0.91, 0.34, 0.26, 0.31, 0.37, 0.45, 0.1, 0.51, 0.22, 0.08, 0.57, 0.39, 0.92, 0.48, 0.99, 0.45, 0.6, 0.94, 0.65, 0.28, 0.51, 0.86, 0.25, 0.72, 0.89, 0.99, 0.92, 0.01, 0.77, 0.34, 0.22, 0.09, 0.23, 0.77, 0.95, 0.14, 0.88, 0.15, 0.13, 0.77, 0.08, 0.53, 0.63, 0.93, 0.41, 0.84, 0.41, 0.85, 0.43, 0.2, 0.09, 0.58, 0.39, 0.95, 0.22, 0.92, 0.5, 0.93, 0.95, 0.69, 0.29, 0.41, 0.79, 0.16, 0.05, 0.66, 0.53, 0.21, 0.93, 0.68, 0.14, 0.96, 0.65, 0.99, 0.1, 0.14, 0.32, 0.82, 0.96, 0.72, 0.7, 0.9, 0.51, 0.15, 0.7, 0.69, 0.93, 0.36, 0.93, 0.93, 0.19, 0.69, 0.81, 0.05, 0.5, 0.34, 0.8, 0.43, 0.72, 0.29, 0.16, 0.63, 0.64, 0.56, 0.91, 0.47, 0.02, 0.73, 0.88, 0.74, 0.7, 0.83, 0.97, 0.96, 0.17, 0.59, 0.07, 0.32, 0.43, 0.83, 0.17, 0.16, 0.35, 0.56, 0.24, 0.38, 0.1, 0.68, 0.61, 0.54, 0.1, 0.43, 0.95, 0.96, 0.92, 0.79, 0.26, 0.79, 0.95, 0.37, 0.5, 0.28, 0.74, 0.89, 0.96, 0.4, 0.39, 0.97, 0.17, 0.12, 0.88, 0.03, 0.9, 0.8, 0.74, 0.1, 0.81, 0.79, 0.28, 0.68, 0.36, 0.17, 0.0, 0.01, 0.21, 0.57, 0.85, 0.15, 0.95, 0.34, 0.01, 0.22, 0.44, 0.38, 0.31, 0.93, 0.0, 0.33, 0.08, 0.51, 0.23, 0.37, 0.71, 1.0, 1.0, 0.4, 0.05, 0.11, 0.0, 0.0, 0.64, 0.08, 0.48, 0.34, 0.4, 0.39, 0.68, 0.43, 0.25, 0.9, 0.0, 0.0, 0.0, 0.41, 0.63, 0.74, 0.87, 0.02, 0.63, 0.36, 0.07, 0.41, 0.0, 0.0, 0.0, 0.0, 0.31, 0.87, 0.31, 0.73, 0.55, 0.46, 0.05, 0.93, 0.0, 0.0, 0.0, 0.0, 0.0, 0.58, 0.24, 0.43, 0.89, 0.18, 0.25, 0.44, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.49, 0.71, 0.14, 0.24, 0.61, 0.61, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.71, 0.29, 0.34, 0.76, 0.21, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.96, 0.9, 0.84, 0.35, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.68, 0.98, 0.12, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.03, 0.17, 0.52, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.3, 0.3, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.3, 0.0, 0.27, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.3, 0.27, 0.0};

const int cupct2_g[CUPCT2_ADJSIZE] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0, 0, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 1, 0, 1, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0};