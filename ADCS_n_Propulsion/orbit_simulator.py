import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
import matplotlib.pyplot as plt

def drag_calc(v_mag, r_pos):
    a_frontal = 1
    c_d = 2.2
    # get drag on r_pos
    rho = 1.38E-11
    return 0.5*rho**c_d*a_frontal*v_mag**2

def dy_vector(t, Y, mu, g0, Isp, T, t_burn, r_goal, burn_r, r_body):
    # Y state vector
    # first 2 x and y location
    # third and fourth are dx and dy speed
    # last is mass
    # Dy is therefore
    # first 2 are dx and dy speed
    # third and fourth are d**2x d**2y acceleration
    # last is dmass
    dy = np.zeros(5)
    p_vec = Y[:2]  # position vector
    r_pos = np.linalg.norm(p_vec)  # radius of position

    v_vec = Y[2:4]  # speed vector
    v_mag = np.linalg.norm(v_vec)  # magnitude of the total speed

    m = Y[4]  # current mass

    dy[:2] = v_vec   # speed vector
    # acceleration f = ma, m1*m2*G/r^2 = m*a, mu*m /r^2 = ma (mu = G*m1), mu /r^2 = a

    dy[2:4] = -mu * Y[:2] / (r_pos**3)
    f_aero = drag_calc(v_mag, r_pos)
    dy[2:4] += v_vec / v_mag * f_aero / m

    if t < t_burn and r_pos > burn_r:
        dy[2:4] += v_vec/v_mag * T / m
        dy[4] = - np.abs(T)/(Isp * g0)  # abs to get reduction always even if thrusts decelerates

    return dy

def reached_orbit(t, Y, mu, g0, Isp, T, t_burn, r_goal, burn_r, r_body):
    """Determine if the spacecraft reaches the destination radius.

    Returns the difference between the current orbital radius and the
    destination radius.
    """
    r_vec = Y[0:2]
    r = np.linalg.norm(r_vec)
    # event only does something if sign changes, if r_goal is infinite never changes
    return r - r_goal

reached_orbit.terminal = True

def run_orbital_maneuver(mu=398600.435507*1e9, g0=9.80655, r_body=6_378_136.6, Isp=300, r_sat=300_000, m_sat=50,
                         T_sat=0, Y_0=False, t_end=365*24*60*60, t_burn=np.inf, r_goal=np.inf, burn_r=0):
    r_start = r_body + r_sat  #
    v_orbit = (mu / r_start) ** 0.5  # m/s, verified
    # t_orbit = 2 * np.pi * r_start / v_orbit

    t_eval = np.linspace(0, t_end, int(1E6))
    if not Y_0:
        Y_0 = np.array([0, r_start, v_orbit, 0, m_sat])

    sol = solve_ivp(
        dy_vector,
        t_span=(0, t_end),
        y0=Y_0,
        t_eval=t_eval,
        args=(mu, g0, Isp, T_sat, t_burn, r_goal+r_body, burn_r+r_body, r_body),
        events=reached_orbit,
        rtol=1E-12,
        atol=1E-15,
        method="DOP853",  # fortran can have quick errors already in *10-2
    )

    return Y_0, sol.y.T[-1], sol

def unit_tests():
    def unit_test1():
        # test if a non-disturbed input has the same output if one orbital period is elapsed

        mu = 3.986004418e14  # earth
        r_earth = 6_378_136.6
        r_sat = 500_000
        r_start = r_earth + r_sat  #
        v_orbit = (mu / r_start) ** 0.5  # m/s, verified
        t_orbit = 2 * np.pi * r_start / v_orbit

        Y_0, Y_end, Y_sol = run_orbital_maneuver(t_end=t_orbit, mu=mu, r_body=r_earth, r_sat=r_sat)
        assert(np.all(np.isclose(Y_0, Y_end, atol=1E-4)))


    def unit_test2():
        # test if a non-disturbed input has an opposite sign radius and velocity input for half an orbital period

        mu = 3.986004418e14  # earth
        r_earth = 6_378_136.6
        r_sat = 500_000
        r_start = r_earth + r_sat  #
        v_orbit = (mu / r_start) ** 0.5  # m/s, verified
        t_orbit = 2 * np.pi * r_start / v_orbit

        Y_0, Y_end, sol = run_orbital_maneuver(t_end=t_orbit/2, mu=mu, r_body=r_earth, r_sat=r_sat)
        Y_test = np.array([-Y_end[0], -Y_end[1], -Y_end[2], -Y_end[3], Y_end[4]])
        assert(np.all(np.isclose(Y_0, Y_test, atol=1E-4)))


    def unit_test3():
        # test if a non-disturbed input shifts r_pos to x and x_vel to -y_vel for 1/4 of an orbital period

        mu = 3.986004418e14  # earth
        r_earth = 6_378_136.6
        r_sat = 500_000
        r_start = r_earth + r_sat  #
        v_orbit = (mu / r_start) ** 0.5  # m/s, verified
        t_orbit = 2 * np.pi * r_start / v_orbit

        Y_0, Y_end, sol = run_orbital_maneuver(t_end=t_orbit / 4, mu=mu, r_body=r_earth, r_sat=r_sat)
        Y_test = [Y_0[1], Y_0[0], Y_0[3], -Y_0[2], Y_0[4]]
        assert(np.all(np.isclose(Y_test, Y_end, atol=1E-4)))

    def unit_test4():
        # test if non disturbed input for circular input doesn't decrease in height
        mu = 3.986004418e14  # earth
        r_earth = 6_378_136.6
        r_sat = 300_000
        r_start = r_earth + r_sat  #
        v_orbit = (mu / r_start) ** 0.5  # m/s, verified
        t_orbit = 2 * np.pi * r_start / v_orbit

        Y_0, Y_end, sol = run_orbital_maneuver(t_end=t_orbit)
        min_distance = r_sat
        for Y in sol.y.T:
            r_distance = (Y[0]**2 + Y[1]**2)**0.5 - r_earth
            if min_distance > r_distance:
                min_distance = r_distance

        assert(np.isclose(min_distance, r_sat, atol=1E-4))  # 10 meter difference allowed

        # varying mass
        Y_0, Y_end, sol = run_orbital_maneuver(t_end=t_orbit, m_sat=10000)
        min_distance = r_sat
        for Y in sol.y.T:
            r_distance = (Y[0]**2 + Y[1]**2)**0.5 - r_earth
            if min_distance > r_distance:
                min_distance = r_distance

        assert(np.isclose(min_distance, r_sat, atol=1))  # 10 meter difference allowed

    def unit_test5():
        # Based on a calculated hohmann transfer do an impulsive burn calculation
        r_sat = 500_000
        r_goal = 100_000
        r_earth = 6_378_136.6
        mu = 3.986004418e14
        v_orbit_start = (mu / (r_sat+r_earth)) ** 0.5  # m/s, verified
        v_orbit_goal = (mu*((2/(r_earth+r_sat))-(2/(2*r_earth+r_goal+r_sat))))**0.5

        T_hohmann_sat = -150000000
        m_sat = 500
        dV = np.abs(v_orbit_start-v_orbit_goal)
        g0 = 9.80665
        Isp = 300
        t_burn_hohmann = np.abs((1-(1/np.exp(dV/(g0*Isp))))*g0*Isp*m_sat/T_hohmann_sat)
        # regular hohmann transfer 500 to 100 km, satellite 500 kg
        Y_0, Y_end, sol = run_orbital_maneuver(t_burn=t_burn_hohmann, r_sat=r_sat, Isp=Isp, g0=g0, mu=mu,
                                               T_sat=T_hohmann_sat, m_sat=m_sat, t_end=60*60*5, r_goal=150_000)
        distance = []
        min_distance = r_sat
        for Y in sol.y.T: # this doesn't actually fully work? Its weird but when inspected this doesn't actually get the minimum number
            r_distance = (Y[0]**2 + Y[1]**2)**0.5 - r_earth
            distance.append(r_distance/1000)
            if min_distance > r_distance:
                min_distance = r_distance

        print(min_distance)
        plt.plot(sol.t, distance)
        plt.show()
        print(r_goal)
        assert (np.isclose(r_goal, min_distance, atol=1E-4))

    unit_test1()
    unit_test2()
    unit_test3()
    unit_test4()
    #unit_test5()



if __name__ == "__main__":
    unit_tests()

    def sat_run():
        # test if a non-disturbed input shifts r_pos to x and x_vel to -y_vel for 1/4 of an orbital period
        r_goal = (6371+200)*1e3
        Y_0, Y_end, sol = run_orbital_maneuver(r_goal=r_goal)

