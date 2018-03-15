"""
Test distance formulas
"""
import math
import unittest

from geopy.point import Point
from geopy.distance import (Distance,
                            GreatCircleDistance,
                            VincentyDistance,
                            EARTH_RADIUS,
                            ELLIPSOIDS)


EARTH_CIRCUMFERENCE = 2 * math.pi * EARTH_RADIUS
NORTH_POLE = Point(90, 0)
SOUTH_POLE = Point(-90, 0)
FIJI = Point(-16.1333333, 180.0) # Vunikondi, Fiji


class CommonDistanceComputationCases(object):

    cls = None

    def test_zero_measure(self):
        self.cls(
            (40.753152999999998, -73.982275999999999),
            (40.753152999999998, -73.982275999999999)
        )

    def test_should_have_length_when_only_given_length(self):
        distance = 1
        self.assertEqual(self.cls(distance).kilometers, distance)

    def test_should_have_zero_distance_for_coincident_points(self):
        self.assertEqual(self.cls((0, 0), (0, 0)).kilometers, 0)

    def test_should_have_nonzero_distance_for_distinct_points(self):
        self.assertTrue(self.cls((0, 0), (0, 1)).kilometers > 0)

    def test_max_longitude(self):
        distance = self.cls(kilometers=1.0)
        destination = distance.destination(FIJI, 45)
        self.assertAlmostEqual(destination.longitude, -179.99338, 4)

    def test_should_compute_distance_for_trip_between_poles(self):
        distance = self.cls(SOUTH_POLE, NORTH_POLE)
        expected_distance = EARTH_CIRCUMFERENCE / 2
        self.assertAlmostEqual(distance.kilometers, expected_distance, -2)

    def test_should_compute_destination_for_trip_between_poles(self):
        distance = self.cls(EARTH_CIRCUMFERENCE / 2)
        destination = distance.destination(NORTH_POLE, 0)
        self.assertAlmostEqual(destination.latitude, -90, 0)
        self.assertAlmostEqual(destination.longitude, 0)

    def test_should_recognize_equivalence_of_pos_and_neg_180_longitude(self):
        distance = self.cls((0, 180), (0, -180)).kilometers
        self.assertAlmostEqual(distance, 0)


class CommonMathematicalOperatorCases(object):

    cls = None

    def test_should_be_able_to_add_distances(self):
        added = self.cls(1.0) + self.cls(1.0)
        self.assertAlmostEqual(added.kilometers, 2.0)

    def test_should_not_allow_adding_with_objects_that_arent_distances(self):
        with self.assertRaises(TypeError):
            self.cls(1.0) + 5

    def test_should_be_able_to_negate_distances(self):
        distance = self.cls(1.0)
        self.assertAlmostEqual(-(distance.kilometers),
                               (-distance).kilometers)

    def test_should_be_able_to_subtract_distances(self):
        subtracted = self.cls(2.0) - self.cls(1.0)
        self.assertAlmostEqual(subtracted.kilometers, 1)

    def test_should_be_able_to_multiply_distances_by_floats(self):
        self.assertAlmostEqual((self.cls(2.0) * 2.0).kilometers, 4.0)

    def test_should_not_be_able_to_multiply_distances_by_distances(self):
        with self.assertRaises(TypeError):
            self.cls(1.0) * self.cls(2.0)

    def test_should_be_able_to_divide_distances_by_distances(self):
        ratio = self.cls(4.0) / self.cls(2.0)
        self.assertAlmostEqual(ratio, 2.0)

    def test_should_be_able_to_divide_distances_by_floats(self):
        divided_distance = self.cls(4.0) / 2.0
        self.assertAlmostEqual(divided_distance.kilometers, 2.0)

    def test_should_be_able_to_take_absolute_value_of_distances(self):
        self.assertAlmostEqual(abs(self.cls(-1.0)).kilometers, 1.0)

    def test_should_be_true_in_boolean_context_when_nonzero_length(self):
        self.assertTrue(self.cls(1.0))

    def test_should_be_false_in_boolean_context_when_zero_length(self):
        self.assertFalse(self.cls(0))

    def test_should_get_consistent_results_for_distance_calculations(self):
        distance1, distance2 = [self.cls((0, 0), (0, 1))
                                for _ in range(2)]
        self.assertEqual(distance1.kilometers, distance2.kilometers)


class CommonConversionCases(object):

    cls = None

    def test_should_convert_to_kilometers(self):
        self.assertEqual(self.cls(1.0).kilometers, 1.0)

    def test_should_convert_to_kilometers_with_abbreviation(self):
        self.assertEqual(self.cls(1.0).km, 1.0)

    def test_should_convert_to_meters(self):
        self.assertEqual(self.cls(1.0).meters, 1000.0)

    def test_should_convert_to_meters_with_abbreviation(self):
        self.assertEqual(self.cls(1.0).m, 1000.0)

    def test_should_convert_to_miles(self):
        self.assertAlmostEqual(self.cls(1.0).miles, 0.621371192)

    def test_should_convert_to_miles_with_abbreviation(self):
        self.assertAlmostEqual(self.cls(1.0).mi, 0.621371192)

    def test_should_convert_to_feet(self):
        self.assertAlmostEqual(self.cls(1.0).feet, 3280.8399, 4)

    def test_should_convert_to_feet_with_abbreviation(self):
        self.assertAlmostEqual(self.cls(1.0).ft, 3280.8399, 4)

    def test_should_convert_to_nautical_miles(self):
        self.assertAlmostEqual(self.cls(1.0).nautical, 0.539956803)

    def test_should_convert_to_nautical_miles_with_abbrevation(self):
        self.assertAlmostEqual(self.cls(1.0).nm, 0.539956803)

    def test_should_convert_from_meters(self):
        self.assertEqual(self.cls(meters=1.0).km, 0.001)

    def test_should_convert_from_feet(self):
        self.assertAlmostEqual(self.cls(feet=1.0).km, 0.0003048, 8)

    def test_should_convert_from_miles(self):
        self.assertAlmostEqual(self.cls(miles=1.0).km, 1.6093440)

    def test_should_convert_from_nautical_miles(self):
        self.assertAlmostEqual(self.cls(nautical=1.0).km, 1.8520000)


class CommonComparisonCases(object):

    cls = None

    def test_should_support_comparison_with_distance(self):
        self.assertTrue(self.cls(1) <= self.cls(1))
        self.assertTrue(self.cls(1) >= self.cls(1))
        self.assertTrue(self.cls(1) < self.cls(2))
        self.assertTrue(self.cls(2) > self.cls(1))
        self.assertTrue(self.cls(1) == self.cls(1))
        self.assertTrue(self.cls(1) != self.cls(2))

        self.assertFalse(self.cls(2) <= self.cls(1))
        self.assertFalse(self.cls(1) >= self.cls(2))
        self.assertFalse(self.cls(2) < self.cls(1))
        self.assertFalse(self.cls(1) > self.cls(2))
        self.assertFalse(self.cls(1) == self.cls(2))
        self.assertFalse(self.cls(1) != self.cls(1))

    def test_should_support_comparison_with_number(self):
        self.assertTrue(1 <= self.cls(1))
        self.assertTrue(1 >= self.cls(1))
        self.assertTrue(1 < self.cls(2))
        self.assertTrue(2 > self.cls(1))
        self.assertTrue(1 == self.cls(1))
        self.assertTrue(1 != self.cls(2))

        self.assertFalse(2 <= self.cls(1))
        self.assertFalse(1 >= self.cls(2))
        self.assertFalse(2 < self.cls(1))
        self.assertFalse(1 > self.cls(2))
        self.assertFalse(1 == self.cls(2))
        self.assertFalse(1 != self.cls(1))


class CommonDistanceCases(CommonDistanceComputationCases,
                          CommonMathematicalOperatorCases,
                          CommonConversionCases,
                          CommonComparisonCases):
    pass


class TestWhenInstantiatingBaseDistanceClass(unittest.TestCase):
    def test_should_not_be_able_to_give_multiple_points(self):
        with self.assertRaises(NotImplementedError):
            Distance(1, 2, 3, 4)


class TestWhenComputingGreatCircleDistance(CommonDistanceCases,
                                           unittest.TestCase):
    cls = GreatCircleDistance

    def test_should_compute_distance_for_half_trip_around_equator(self):
        distance_around_earth = self.cls((0, 0), (0, 180)).kilometers
        self.assertEqual(distance_around_earth, EARTH_CIRCUMFERENCE / 2)

    def test_should_compute_destination_for_half_trip_around_equator(self):
        distance = self.cls(EARTH_CIRCUMFERENCE / 2)
        destination = distance.destination((0, 0), 0)
        self.assertAlmostEqual(destination.latitude, 0)
        self.assertAlmostEqual(destination.longitude, 180)


class TestWhenComputingVincentyDistance(CommonDistanceCases,
                                        unittest.TestCase):

    cls = VincentyDistance

    def setUp(self):
        self.original_ellipsoid = self.cls.ELLIPSOID

    def tearDown(self):
        self.cls.ELLIPSOID = self.original_ellipsoid

    def test_should_not_converge_for_half_trip_around_equator(self):
        with self.assertRaises(ValueError):
            self.cls((0, 0), (0, 180))

    def test_should_compute_destination_for_half_trip_around_equator(self):
        distance = self.cls(EARTH_CIRCUMFERENCE / 2)
        destination = distance.destination((0, 0), 0)
        self.assertAlmostEqual(destination.latitude, 0, 0)
        self.assertAlmostEqual(destination.longitude, -180, 0)

    def test_should_compute_same_destination_as_other_libraries(self):
        distance = self.cls(54.972271)
        destination = distance.destination((-37.95103, 144.42487), 306.86816)
        self.assertAlmostEqual(destination.latitude, -37.6528177174, 10)
        self.assertAlmostEqual(destination.longitude, 143.9264976682, 10)

    def test_should_get_distinct_results_for_different_ellipsoids(self):
        results = [
            self.cls((0, 0), (0, 1), ellipsoid=ELLIPSOIDS[ellipsoid_name])
            for ellipsoid_name in ELLIPSOIDS.keys()
        ]

        self.assertFalse(any(results[x].kilometers == results[y].kilometers
                             for x in range(len(results))
                             for y in range(len(results))
                             if x != y))
