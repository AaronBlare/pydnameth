import unittest
from tests.definitions import ROOT_DIR
from pydnameth import Data
from pydnameth import Experiment
from pydnameth import Annotations
from pydnameth import Observables
from pydnameth import Cells
from pydnameth import Attributes
from pydnameth import Config
from pydnameth.infrastucture.load.attributes import load_attributes_dict
from pydnameth.infrastucture.load.attributes import load_cells_dict
from pydnameth.config.attributes.subset import pass_indexes
from pydnameth.config.attributes.subset import get_indexes
from pydnameth.config.attributes.subset import subset_attributes
from pydnameth.config.attributes.subset import subset_cells
from pydnameth.config.common import CommonTypes

class TestLoadAnnotations(unittest.TestCase):

    def setUp(self):

        data = Data(
            name='cpg_beta',
            path=ROOT_DIR,
            base='fixtures'
        )

        experiment = Experiment(
            type=None,
            task=None,
            method=None,
            params=None
        )

        annotations = Annotations(
            name='annotations',
            exclude='none',
            cross_reactive='ex',
            snp='ex',
            chr='NS',
            gene_region='yes',
            geo='any',
            probe_class='any'
        )

        observables = Observables(
            name='observables',
            types={}
        )

        cells = Cells(
            name='cells',
            types='any'
        )

        attributes = Attributes(
            target='age',
            observables=observables,
            cells=cells
        )

        self.config = Config(
            data=data,
            experiment=experiment,
            annotations=annotations,
            attributes=attributes,
            is_run=True,
            is_root=True
        )

    def test_pass_indexes_num_elems(self):
        self.config.attributes.observables.types = {'gender': 'any'}
        self.config.attributes_dict = load_attributes_dict(self.config)
        indexes = pass_indexes(self.config, 'gender', 'any', 'any')
        self.assertEqual(len(indexes), 729)

    def test_pass_indexes_num_f(self):
        self.config.attributes.observables.types = {'gender': 'F'}
        self.config.attributes_dict = load_attributes_dict(self.config)
        indexes = pass_indexes(self.config, 'gender', 'F', 'any')
        self.assertEqual(len(indexes), 388)

    def test_pass_indexes_num_m(self):
        self.config.attributes.observables.types = {'gender': 'M'}
        self.config.attributes_dict = load_attributes_dict(self.config)
        indexes = pass_indexes(self.config, 'gender', 'M', 'any')
        self.assertEqual(len(indexes), 341)

    def test_get_indexes_num_elems(self):
        self.config.attributes.observables.types = {'gender': 'any'}
        self.config.attributes_dict = load_attributes_dict(self.config)
        indexes = get_indexes(self.config)
        self.assertEqual(len(indexes), 729)

    def test_get_indexes_num_f(self):
        self.config.attributes.observables.types = {'gender': 'F'}
        self.config.attributes_dict = load_attributes_dict(self.config)
        indexes = get_indexes(self.config)
        self.assertEqual(len(indexes), 388)

    def test_get_indexes_num_m(self):
        self.config.attributes.observables.types = {'gender': 'M'}
        self.config.attributes_dict = load_attributes_dict(self.config)
        indexes = get_indexes(self.config)
        self.assertEqual(len(indexes), 341)

    def test_get_indexes_num_f_m(self):
        self.config.attributes.observables.types = {'gender': ['F', 'M']}
        self.config.attributes_dict = load_attributes_dict(self.config)
        indexes = get_indexes(self.config)
        self.assertEqual(len(indexes), 729)

    def test_subset_attributes(self):
        self.config.attributes_dict = load_attributes_dict(self.config)
        self.config.attributes_indexes = list(range(5))
        subset_attributes(self.config)
        self.assertEqual(self.config.attributes_dict['gender'], ['M']*5)

    def test_subset_cells(self):
        self.config.cells_dict = load_cells_dict(self.config)
        self.config.attributes_indexes = list(range(5))
        subset_cells(self.config)
        self.assertEqual(self.config.cells_dict['CD8T'], [0, 0, 0.006011666, 0, 0])

    def test_incorrect_variable(self):
        self.config.attributes_dict = load_attributes_dict(self.config)
        self.assertRaises(ValueError, pass_indexes, self.config, 'age', 100, CommonTypes.any.value)

    def test_len_list_pass_indexes(self):
        self.config.attributes_dict = load_attributes_dict(self.config)
        self.assertEqual(len(pass_indexes(self.config, 'age', 18, CommonTypes.any.value)), 19)

    def test_incorrect_observables(self):
        self.config.attributes_dict = load_attributes_dict(self.config)
        self.config.attributes.observables.types = {'some_obs': ''}
        self.assertRaises(ValueError, get_indexes, self.config)

    def test_get_indexes_age_and_m(self):
        self.config.attributes_dict = load_attributes_dict(self.config)
        self.config.attributes.observables.types = {'age': [20, 21, 22], 'gender': 'M'}
        self.assertEqual(len(get_indexes(self.config)), 14)


if __name__ == '__main__':
    unittest.main()
