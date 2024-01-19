#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  7 16:38:50 2021

This file imports and checks the Particle3D class by performing some basic
tests. Amongst the things tested are:
    * basic instantiation
    * instance methods
    * creating instances from a data file
    * static methods

@author: miguel and joe
"""
import numpy as np
import pytest
import types


try:
    import particle3D
except SyntaxError as err:
    print("----------------------------------------------------------------------")
    print("There is a syntax error (something that is not valid python code)")
    print("in your exercise1.py.  The message below should help you find it.")
    print(f"Note that the issue may be just above the line shown ({err.lineno}).")
    print("----------------------------------------------------------------------")
    raise
except ImportError:
    raise ImportError(
        "You have not created particle3D.py in the same directory as this file, "
        "or there was a syntax error in it."
    )


def test_exists():
    assert hasattr(particle3D, "Particle3D"), "There is no Particle3D class in particle3D.py"
    assert type(particle3D.Particle3D) == type, "Particle3D in particle3D.py is not a class - either you overwrote it or something else"

def test_init():
    m0 = 1.0
    x0 = [1.0, 2.0, 3.0]
    v0 = np.array([1.0, -1.0, 0.0])
    p = particle3D.Particle3D("new", m0, x0, v0)
    assert hasattr(p, "label"), "Newly-created particles have no attribute 'label'"
    assert hasattr(p, "mass"), "Newly-created particles have no attribute 'mass'"
    assert hasattr(p, "position"), "Newly-created particles have no attribute 'position'"
    assert hasattr(p, "velocity"), "Newly-created particles have no attribute 'velocity'"

    assert p.label == "new", "Particle label not set correctly"
    assert np.isclose(p.mass, 1.0,), "Particle mass not set correctly"
    assert isinstance(p.position, np.ndarray), "Particle positon not converted into an array from a list"
    assert np.allclose(p.position, [1.0, 2.0, 3.0]), "Particle positon not set correctly"

    # Not tested: velocity

def test_kinetic_energy():
    m0 = 4.0
    x0 = np.array([1.0, 2.0, 3.0])
    v0 = np.array([-2.0, 4.0, 0])
    p = particle3D.Particle3D("new", m0, x0, v0)

    assert hasattr(p, "kinetic_energy"), "Class has no kinetic_energy method"
    assert type(p.kinetic_energy) == types.MethodType, "kinetic_energy is not a proper method"
    k = p.kinetic_energy()
    assert k is not None, "Kinetic energy returned no value"
    assert np.isclose(k, 40.0), "Kinetic energy gives wrong value"


def test_position_update():
    m0 = 1.0
    x0 = np.array([1.0, 2.0, 3.0])
    v0 = np.array([1.0, -1.0, 0])
    p = particle3D.Particle3D("new", m0, x0, v0)

    dt = 0.125
    expected_position_change = dt * np.array([1, -1, 0])

    start_position = p.position.copy()
    p.update_position_1st(dt)
    position_change = p.position - start_position

    assert np.allclose(position_change, expected_position_change), "update_position_1st did not move particle by correct amount"


def test_velocity_update():
    m0 = 2.0
    x0 = np.array([1.0, 2.0, 3.0])
    v0 = np.array([1.0, -1.0, 0])
    p = particle3D.Particle3D("new", m0, x0, v0)

    dt = 0.125
    force = np.array([-1.0, 0.0, 2.0])
    expected_velocity_change = dt * force / m0

    start_velocity = p.velocity.copy()
    p.update_velocity(dt, force)
    velocity_change = p.velocity - start_velocity

    assert np.allclose(velocity_change, expected_velocity_change), "update_velocity did not accelerate particle by correct amount"


def test_read_line():
    p_list = []
    for line in open("p3d_test_data.txt"):
        line = line.strip()
        if line.startswith('#') or line == "":
            continue
        try:
            p = particle3D.Particle3D.read_line(line)
        except Exception as error:
            raise AssertionError(f"read_line function crashed when reading the line '{line}' with this message: {error}")
        p_list.append(p)

    assert p_list[0] is not None, "read_file returned no value"
    assert isinstance(p_list[0], particle3D.Particle3D), "read_file did not return a Particle3D instance"
    assert p_list[0].label == "Alice", "First particle label not correct"
    assert p_list[1].label == "Bob", "Second particle label not correct"
    assert p_list[2].label == "Charlie", "Third particle label not correct"

    assert not isinstance(p_list[0].mass, str), "You did not convert mass to float in read_line - it is a string"

    assert np.isclose(p_list[0].mass, 1.0), "First particle mass not correct in read_line"
    assert np.isclose(p_list[1].mass, 0.5), "Second particle mass not correct in read_line"
    assert np.isclose(p_list[2].mass, -1.0), "Third particle mass not correct in read_line"

    assert not isinstance(p_list[0].position[0], str), "You did not convert position to a float array in read_line - it is a string"
    assert not isinstance(p_list[0].velocity[0], str), "You did not convert velocity to a float array in read_line - it is a string"

    assert np.allclose(p_list[0].position, [0, 0, 0]), "First particle position not correct"
    assert np.allclose(p_list[1].position, [0.5, -0.5, 0]), "Second particle position not correct"
    assert np.allclose(p_list[2].position, [0, 0, 1.0]), "Third particle position not correct"

    # TODO Check velocity


def test_total_kinetic():
    p_list = [
        particle3D.Particle3D("A", 1.0, np.array([0.0, 0.0, 0.0]), np.array([0.0, 0.0, 1.0])),
        particle3D.Particle3D("B", 0.5, np.array([0.0, 0.0, 0.0]), np.array([0.5, 0.5, 0.0])),
        particle3D.Particle3D("C", 2.0, np.array([0.0, 0.0, 0.0]), np.array([-1.0, 1.0, 0.0])),
    ]

    total_ke = particle3D.Particle3D.total_kinetic_energy(p_list)
    assert total_ke is not None, "Total kinetic energy did not return a value"
    assert np.isclose(total_ke, 2.625), "Total kinetic energy incorrect"


def test_str():
    p = particle3D.Particle3D("Alice", 1.0, np.array([1.0, 2.0, 3.0]), np.array([-1.0, -2.0, -3.0]))

    s = str(p)

    words = s.split()
    assert "\n" not in s, "There should be no new line in str method"
    assert len(words) == 4, "str method should return string with 4 words in it"
    assert words[0] == "Alice", "str method label not correct"
    assert eval(words[1]) == 1.0, "str method x not correct"
    assert eval(words[2]) == 2.0, "str method y not correct"
    assert eval(words[3]) == 3.0, "str method z not correct"


if __name__ == "__main__":
    status = pytest.main(["-s", __file__])
    if status == 0:
        print(
            "All of your tests pass, but this does not mean automatic full marks as we don't test everything."
        )
    else:
        print(
            "There are one or more things wrong with your code. See the messages above, and the details above that."
        )
