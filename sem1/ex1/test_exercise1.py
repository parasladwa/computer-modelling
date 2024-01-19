"""
You don't have to understand this file unless you want to.

This module contains what are called "Unit tests" - that is,
tests of individual functions that you have written to make
sure they give the right answer.

Real code should contain tests! In industry it is almost universal;
many companies or coders write the test code first and then fill in
the code that makes the test pass, as you are doing here.

The pytest program looks for functions that start or end with the
word "test", and runs each of them in a separate python program.

These test are not written quite like tests in real code, because
I have designed them to give you more specific guidance when they
fail. That can make them a bit harder to read, but what they print
out will be more helpful.

They are also not quite complete - see the bottom of test_task2 for
details.
"""
print("This is version 1 of test_exercise1.py")
import pytest
import numpy as np

try:
    import exercise1
except SyntaxError as err:
    print("----------------------------------------------------------------------")
    print("There is a syntax error (something that is not valid python code)")
    print("in your exercise1.py.  The message below should help you find it.")
    print(f"Note that the issue may be just above the lineshown ({err.lineno}).")
    print("----------------------------------------------------------------------")
    raise
except ImportError:
    raise ImportError(
        "You have not created exercise1.py in the same directory as this file, "
        "or there was a syntax error in it."
    )


def test_task1():
    # check that you did not rename your function
    # We use the "assert" statement in tests to indicate that if everything is
    # working the first expression should be true. If it is not true then the
    # error message in the second part is printed instead.
    # The python function "hasattr" checks whether (in this case), exercise1
    # has something (a function) in it called "task2".
    assert hasattr(exercise1, "task1"), "task1 has been renamed or not present!"

    # for b = True
    # we are doing 3^1 + 4^1 + 5^1, which is 3 + 4 + 5 = 12
    x = [3.0, 4.0, 5.0]
    n = 1
    b = True
    expected = 12.0
    try:
        result = exercise1.task1(x, n, b)
    except:
        raise RuntimeError(f"task1 crashed when called with arguments ({x}, {n}, {b})")

    # This assertion will cause an error if the first term is not true, which will
    # cause pytest to mark this test as failed.

    # We use np.isclose because floating point arithmetic is not exact, because values
    # are stored to a finite number of decimal places.
    assert np.isclose(
        result, expected
    ), f"task1({x}, {n}, {b}) should be 50 but is {result}"

    # Now we test for b = False.
    # we are doing 3^2 + 4^2 + 5^2, which is 9 + 16 + 25 = 50

    b = False
    expected = 50.0
    try:
        result = exercise1.task1(x, n, b)
    except:
        raise RuntimeError(f"task1 crashed when called with arguments ({x}, {n}, {b})")

    assert np.isclose(
        result, expected
    ), f"task1({x}, {n}, {b}) should be 50 but is {result}"


def test_task2():
    # Test that you have written your function
    assert hasattr(exercise1, "task2"), "task2 not written yet"

    # Try running your function with two arguments and see if it works at all
    try:
        result = exercise1.task2(2, 3)
    except:
        raise RuntimeError("task2 crashed when called with arguments (2, 3)")

    # Check that you did not forget the "return ..." at the end of your function
    assert result is not None, "task2 did not return anything at all"

    # Check that you returned multiple things, not just one
    try:
        number_of_results = len(result)
    except TypeError:
        raise ValueError("task2 only returned a single value not three things")

    # Check that you returned three things
    assert number_of_results == 3, "task2 did not return three things"

    # Check that all three of the things returned are arrays
    for i in range(3):
        assert (
            type(result[i]) == np.ndarray
        ), f"task2 did not return an array for its output {i}"

    # Check that the first result is an array of all zeros
    a = result[0]
    assert a.shape == (2,), "task2 first result is the wrong shape"
    assert (a == 0).all(), "task2 first result should be all zeros but is not"

    assert exercise1.task2.__doc__ is not None, "task2 has no docstring"

    # We have not checked the other results here.

    # If you're reading this then you are exploring the code, which is exactly the
    # thing you should be doing. To make sure your code is working you could
    # write additional tests in this function in the places
    # where I have written TODO below.

    # TODO
    # Add a test here to check that result[1] and result[2] are as expected.


def test_task3a():
    # Test that you have written your function
    assert hasattr(exercise1, "task3a"), "task3a not written yet"

    # Try running your function with two arguments and see if it works at all
    a = np.array([1, 2, 3])
    b = np.array([4, 5, 6])
    t = 2
    try:
        result = exercise1.task3a(a, b, t)
    except:
        raise RuntimeError(f"task3a crashed when called with arguments ({a}, {b}, {t})")

    # Check that you did not forget the "return ..." at the end of your function
    assert result is not None, "task3a did not return anything at all"

    # Check that you return an array not a list
    assert type(result) == np.ndarray, "task3a did not return an array"

    # Check that your return value is the same shape as a and b
    assert result.shape == (3,), "task3a result is the wrong shape"

    # Check the result is correct. It is safe to do this with == because the
    # values are all integers, which can be compared exactly, unlike floats.
    assert np.all(result == np.array([20, 28, 36]))


def test_task3b():
    assert hasattr(exercise1, "task3b"), "task3b not written yet"

    # Try running your function with two arguments and see if it works at all
    a = np.array([1.0, 2.0, 3.0])
    b = np.array([2.0, 3.0, 4.0])

    try:
        result = exercise1.task3b(a, b)
    except:
        raise RuntimeError(f"task3b crashed when called with arguments ({a}, {b})")

    assert isinstance(result, float), "task3b did not return a floating point value"

    # Check that the result is sqrt(3), which is the right answer for the vectors above.
    # Note that we can't just use == because these are floating point numbers with only
    # limited decimal places.  We use np.isclose instead.
    assert np.isclose(result, np.sqrt(3))


def test_task4a():
    assert hasattr(exercise1, "task4a"), "task4a not written yet"

    # Try running your function with two arguments and see if it works at all
    v1 = np.array([1.0, 2.0, 3.0])
    v2 = np.array([2.0, 3.0, 4.0])

    try:
        result = exercise1.task4a(v1, v2)
    except:
        raise RuntimeError(f"task4a crashed when called with arguments ({v1}, {v2})")

    # Check that you did not forget the "return ..." at the end of your function
    assert result is not None, "task4a did not return anything at all"

    # Check that you returned multiple things, not just one
    try:
        number_of_results = len(result)
    except TypeError:
        raise ValueError(f"task4a only returned a single value not two vectors")

    # check that you returned two things
    assert (
        number_of_results == 2
    ), "task4a returned the wrong number of things - should be two vectors"

    # check that the types are correct
    assert (
        type(result[0]) == np.ndarray
    ), "task4a did not return an array as its first output"
    assert (
        type(result[1]) == np.ndarray
    ), "task4a did not return an array as its second output"

    # check that the two sides are the same, and that both are the correct answer
    assert np.allclose(
        result[0], result[1]
    ), "task4a returned different vectors but they should be the same"

    # Check that we get the correct values of v1 X v2
    assert np.allclose(result[0], np.array([-1.0, 2.0, -1.0]))


def test_task4b():
    assert hasattr(exercise1, "task4b"), "task4b not written yet"

    # Try running your function with two arguments and see if it works at all
    v1 = np.array([1.0, 2.0, 3.0])
    v2 = np.array([2.0, 3.0, 4.0])
    v3 = np.array([-1.0, 5.0, 2.0])

    try:
        result = exercise1.task4b(v1, v2, v3)
    except:
        raise RuntimeError(
            f"task4b crashed when called with arguments ({v1}, {v2}, {v3})"
        )

    # Check that you did not forget the "return ..." at the end of your function
    assert result is not None, "task4b did not return anything at all"

    # Check that you returned multiple things, not just one
    try:
        number_of_results = len(result)
    except TypeError:
        raise ValueError(f"task4b only returned a single value not two vectors")

    # check that you returned two things
    assert (
        number_of_results == 2
    ), "task4b returned the wrong number of things - should be two vectors"

    # check that the types are correct
    assert (
        type(result[0]) == np.ndarray
    ), "task4b did not return an array as its first output"
    assert (
        type(result[1]) == np.ndarray
    ), "task4b did not return an array as its second output"

    # check that the two sides are the same, and that both are the correct answer
    assert np.allclose(
        result[0], result[1]
    ), "task4b returned different vectors but they should be the same"

    # TODO
    # Add a test here that the checks if the function produces the correct result
    # for the vectors above.


def test_task4c():
    assert hasattr(exercise1, "task4c"), "task4c not written yet"

    # Try running your function with two arguments and see if it works at all
    v1 = np.array([1.0, 2.0, 3.0])
    v2 = np.array([2.0, 3.0, 4.0])
    v3 = np.array([-1.0, 5.0, 2.0])

    try:
        result = exercise1.task4c(v1, v2, v3)
    except:
        raise RuntimeError(
            f"task4c crashed when called with arguments ({v1}, {v2}, {v3})"
        )

    # Check that you did not forget the "return ..." at the end of your function
    assert result is not None, "task4c did not return anything at all"

    # Check that you returned multiple things, not just one
    try:
        number_of_results = len(result)
    except TypeError:
        raise ValueError(f"task4c only returned a single value not two vectors")

    # check that you returned two things
    assert (
        number_of_results == 2
    ), "task4c returned the wrong number of things - should be two vectors"

    # check that the types are correct
    assert (
        type(result[0]) == np.ndarray
    ), "task4c did not return an array as its first output"
    assert (
        type(result[1]) == np.ndarray
    ), "task4c did not return an array as its second output"

    # check that the two sides are the same, and that both are the correct answer
    assert np.allclose(
        result[0], result[1]
    ), "task4c returned different values but they should be the same"

    # TODO
    # Add a test here that the checks if the function produces the correct result
    # for the vectors above.


def test_task5():
    m1 = 2.0
    m2 = 3.0
    x1 = np.array([0.0, 0.0, 0.0])
    x2 = np.array([2.0, 0.0, 0.0])

    assert hasattr(exercise1, "task5"), "task5 not written yet"

    try:
        result = exercise1.task5(x1, x2, m1, m2)
    except:
        raise RuntimeError(
            f"task5 crashed when called with arguments ({x1}, {x2}, {m1}, {m2})"
        )

    assert result is not None, "task5 did not return anything"
    try:
        number_of_results = len(result)
    except TypeError:
        raise ValueError("task5 only returned a single value not two values")

    assert number_of_results == 2, "task5 did not return two values"

    F, phi = result

    # Check correct values are coming back. My x1 and x2 are separated by distance 2,
    # and m1 * m2 = 6, so we expect
    G = 6.6743e-11
    expected_F = -1.5 * G * np.array([1, 0, 0])
    assert np.allclose(F, expected_F), "Force value not as expected"

    # TODO
    # Add a test for phi


def test_task6a():
    assert hasattr(exercise1, "task6a"), "task6a not written yet"
    n = 3
    try:
        M = exercise1.task6a(n)
    except:
        raise RuntimeError(f"task6a crashed when called with argument ({n})")

    assert isinstance(M, np.ndarray), "task6a did not return an array"
    assert M.shape == (n, n), "task6a did not return the correct shaped array (n x n)"

    for i in range(n):
        for j in range(n):
            assert (
                M[i, j] == i + 2 * j
            ), f"task6a gave gave wrong value for element {i},{j} for n = {n}"


def test_task6b():
    assert hasattr(exercise1, "task6b"), "task6b not written yet"
    n = 3
    try:
        y = exercise1.task6b(n)
    except:
        raise RuntimeError(f"task6b crashed when called with argument ({n})")

    assert isinstance(y, np.ndarray), "task6b did not return an array"
    assert y.shape == (
        n,
    ), "task6b did not return the correct shaped array (length n vector)"

    # TODO: Test y values are as expected compared to analytic result




# This special python phrase checks to see if we have run this file or just imported it.
# If we have run the file then the automatic python variable __name__ will be set to the
# value "__main__". Otherwise it will be the name of the imported module. 
if __name__ == "__main__":
    # We use the python tool "pytest" to run all the tests in this file automatically.
    # pytest automatically searches for any functions starting or ending with the word
    # "test" and tries running each in turn. It then reports on all the ones that fail.
    status = pytest.main(["-vs", __file__])

    # Print a little message depending on whether all the tests pass. As in many computer
    # programs, a return value of 0 is used to indicate that everything worked, and any
    # non-zero value means that something went wrong.
    if status == 0:
        print("All of your tests pass, but this does not mean automatic full marks as we don't test everything.")
    else:
        print("There are one or more things wrong with your code. See the messages above, and the details above that.")
