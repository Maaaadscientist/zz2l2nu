import os
import re

import yaml

class SystDatasetSelector:
    """Selector of datasets affected by systematic variations.

    This class provides means to select simulated datasets affected by
    a specific systematic variation or a group of variations.  Datasets
    with real data are always ignored.
    """

    def __init__(self, config_path):
        """Initialize from a configuration file.
        
        The configuration file provides a mapping from labels of
        systematic uncertainties to masks that select datasets affected
        by those uncertainties.
        """

        with open(config_path) as f:
            self.dataset_masks = yaml.safe_load(f)


    def __call__(self, datasets, requested_syst, skip_nominal=False):
        """Select simulated datasets affected by given uncertainty.

        Arguments:
            datasets:  Iterable with datasets.
            requested_syst:  Label of a systematic variation or a group
                             of them.
            skip_nominal:  Do not include the nominal configuration.
                Only applicable if requested_syst is "" or "all".

        Yield value:
            Tuple consisting of a concrete systematic variation and a
            dataset.

        Requested variations can be specified in different ways.  This
        can be a concrete variation (e.g. "jec_up").  If the direction
        posfix is omitted (e.g. "jec") both up and down variations will
        be produced.  If "all" is given, all registered variations and
        also the nominal configuration are produced.  The nominal
        configuration alone can be requested with "".

        This method is a generator.  At each invocation it yields the
        next pair of a concrete systematic variation and a dataset
        affected by it.
        """

        selected_variations = {}

        if not requested_syst:
            # The nominal configuration
            if not skip_nominal:
                selected_variations[''] = '*'
        elif requested_syst == 'all':
            if not skip_nominal:
                selected_variations[''] = '*'

            for syst, masks in self.dataset_masks.items():
                for direction in ['up', 'down']:
                    var_label = '{}_{}'.format(syst, direction)
                    selected_variations[var_label] = masks
        elif requested_syst in self.dataset_masks:
            # A pair of up and down variations has been requested
            syst = requested_syst
            masks = self.dataset_masks[syst]

            for direction in ['up', 'down']:
                selected_variations['{}_{}'.format(syst, direction)] = masks
        else:
            # The only remaining option is a fully specified variation
            syst, direction = self.split_syst_label(requested_syst)

            if not direction or syst not in self.dataset_masks:
                raise RuntimeError(
                    'Group of systematic variations "{}" '
                    'is not recognized.'.format(requested_syst)
                )

            selected_variations[requested_syst] = self.dataset_masks[syst]

        # Systematic variations affect only simulation
        sim_datasets = [d for d in datasets if d.is_sim]

        for variation, masks in selected_variations.items():
            for dataset in sim_datasets:
                selected = False

                for mask in masks:
                    if mask == '*' or mask in dataset.name:
                        selected = True
                        break

                if selected:
                    yield variation, dataset


    @staticmethod
    def split_syst_label(syst_label):
        """Split label of systematic variation into stem and direction.

        The direction postfix is optional.  If missing, return the whole
        syst_label as the stem.
        """

        match = re.match('(.+)_(up|down)', syst_label)

        if match:
            return match.group(1), match.group(2)
        else:
            return syst_label, None

