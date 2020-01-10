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
            self.syst_configs = yaml.safe_load(f)


    def __call__(self, datasets, requested_syst, skip_nominal=False,
                 combine_weights=False):
        """Select simulated datasets affected by given uncertainty.

        Arguments:
            datasets:  Iterable with datasets.
            requested_syst:  Label of a systematic variation or a group
                             of them.
            skip_nominal:  Do not include the nominal configuration.
                Only applicable if requested_syst is "", "all", or
                "weights".
            combine_weights:  Requests that all weight-based variations
                are processed together in one pass.

        Yield value:
            Tuple consisting of a concrete systematic variation and a
            dataset.

        Requested variations can be specified in different ways.  This
        can be a concrete variation (e.g. "jec_up").  If the direction
        posfix is omitted (e.g. "jec") both up and down variations will
        be produced.  If "all" is given, all registered variations and
        also the nominal configuration are produced.  The nominal
        configuration alone can be requested with "".  Weight-only
        variations can be requested with "weights".

        This method is a generator.  At each invocation it yields the
        next pair of a concrete systematic variation and a dataset
        affected by it.
        """

        selected_variations = {}

        if not requested_syst:
            # The nominal configuration
            if not skip_nominal:
                selected_variations[''] = '*'
        elif requested_syst == 'weights':
            if not skip_nominal:
                selected_variations['weights'] = '*'
        elif requested_syst == 'all':
            if not skip_nominal:
                selected_variations['weights' if combine_weights else ''] = '*'

            for syst, config in self.syst_configs.items():
                if combine_weights and config['is_weight']:
                    continue
                for direction in ['up', 'down']:
                    var_label = '{}_{}'.format(syst, direction)
                    selected_variations[var_label] = config['processes']
        elif requested_syst in self.syst_configs:
            # A pair of up and down variations has been requested
            syst = requested_syst
            masks = self.syst_configs[syst]['processes']

            for direction in ['up', 'down']:
                selected_variations['{}_{}'.format(syst, direction)] = masks
        else:
            # The only remaining option is a fully specified variation
            syst, direction = self.split_syst_label(requested_syst)

            if not direction or syst not in self.syst_configs:
                raise RuntimeError(
                    'Group of systematic variations "{}" '
                    'is not recognized.'.format(requested_syst)
                )

            selected_variations[requested_syst] = \
                self.syst_configs[syst]['processes']

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


mpl_style = {
    'figure.figsize': (6.0, 4.8),
    
    'axes.labelsize':              'large',
    'axes.formatter.use_mathtext': True,
    'axes.formatter.limits':       (-2, 4),
    
    'xtick.top':          True,
    'xtick.direction':    'in',
    'xtick.minor.top':    True,
    'xtick.minor.bottom': True,
    
    'ytick.right':        True,
    'ytick.direction':    'in',
    'ytick.minor.left':   True,
    'ytick.minor.right':  True,
    
    'lines.linewidth':   1.,
    'lines.markersize':  3.,
    'errorbar.capsize':  1.
}

