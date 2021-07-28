// def solve_modelC(parameters, t, u):
//     parameters = 10 ** parameters
//     # annotate parameters
//     y0 = parameters[0]
//     l = parameters[1]
//     a = parameters[2]
//     y_model = np.zeros(len(t))
//     # solve model
//     y_model[0] = y0
//     for idx in range(len(t) - 1):
//         dt = t[idx + 1] - t[idx]
//         m = (u[idx + 1] - u[idx]) / dt
//         b = u[idx]
//         c = y_model[idx] - a * b / l + a * m / l ** 2
//         y_model[idx + 1] = a * b / l - a * m / l ** 2 + a * m * dt / l + c * np.exp(-l * dt)
//     return y_model