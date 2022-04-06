"""
test_check.py
"""

from seekr2.modules import check

def test_toy_check_systems_within_Voronoi_cells_elber(toy_elber_model):
    check.check_systems_within_Voronoi_cells(toy_elber_model)