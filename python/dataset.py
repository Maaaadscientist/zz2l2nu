import os
import re


class Dataset:
    """Describes a dataset.

    An object of this class contains the following fields:
        path:    Path to the dataset definition file.
        name:    String identifying the dataset.
        files:   Paths to input ROOT files included in the dataset.
        is_sim:  Boolean indicating whether the dataset is real data or
            simulation.
        parameters:  Mapping with all parameters extracted from the
            definition file.  Integer and floating-point values are
            converted into the corresponding native representations.

    """

    def __init__(self, path):
        self.path = path
        self.name = self.shorten_name(os.path.basename(path))
        self.parameters = {}
        self.is_sim = None
        self.files = []

        self._from_txt(path)


    def save(self, path):
        """Save to a text file."""
        self._save_txt(path)


    @staticmethod
    def shorten_name(long_name):
        """Convert dataset filename to standard short name.

        Adapted from [1].
        [1] https://gitlab.cern.ch/HZZ-IIHE/hzz2l2nu/blob/27d5f5b1c5b4004751f873f9f893366e84a3b8d0/Tools/prepareAllJobs.py#L275-279
        """
        # Strip prefix like "Bonzais-" or "Baobabs-" and pruner posfix
        match = re.match(r'^.+?-(.+)-\w+Pruner.+$', long_name)

        if not match:
            raise RuntimeError(
                'Unexpected format for catalogue name "{}".'.format(long_name)
            )

        name = match.group(1)

        # Remove everything starting from given substrings.  They are
        # not necessarily present, though.
        name = name.split('_TuneCUETP8M1')[0]
        name = name.split('_13TeV')[0]

        return name


    def _from_txt(self, path):
        """Initialize from a file in the old-style catalogue format."""

        blank_regex = re.compile(r'^\s*$')

        # Regular expression that matches a line containing a path to an
        # input file.  The first group captures the path.
        path_regex = re.compile(r'^\s*(\S+).*$')

        # Regular expression that matches a line with a configuration
        # parameter.  The first group captures the name of the
        # parameter, the second group captures its value.
        parameter_regex = re.compile(r'^\s*[#\*]\s*(.+)\s*:\s*(.*)\s*$')

        reading_header = True

        definition_file = open(path)

        for line in definition_file:
            if blank_regex.match(line):
                continue

            if reading_header:
                match = parameter_regex.match(line)

                if match:
                    name = match.group(1)
                    value = match.group(2)

                    if name in self.parameters:
                        print(
                            'Parameter "{}" specified multiple times in data '
                            'set definition file "{}". Overwriting old value '
                            '"{}" with "{}".'.format(
                                name, path, self.paramters[name], value
                            )
                        )

                    self.parameters[name] = value
                else:
                    reading_header = False

            if not reading_header:
                match = path_regex.match(line)

                if match:
                    self.files.append(match.group(1))
                else:
                    print(
                        'In dataset definition file "{}" failed to parse '
                        'line\n{}Skipping it.'.format(path, line)
                    )

        definition_file.close()

        # Try to convert values of parameters to native types
        for name, value in self.parameters.items():
            try:
                self.parameters[name] = int(value)
            except ValueError:
                try:
                    self.parameters[name] = float(value)
                except ValueError:
                    pass


        # Extract important parameters.  They are removed from the
        # common dictionary.
        data_type = self.parameters.pop('data type')

        if data_type == 'data':
            self.is_sim = False
        elif data_type == 'mc':
            self.is_sim = True
        else:
            raise RuntimeError(
                'Illegal value "{}" for parameter "data type" found in data '
                'set defition file "{}".'.format(data_type, path)
            )


    def _save_txt(self, path):
        """Save in the old-style catalogue format.

        The order of parameters in the initial file may not be
        reproduced.
        """

        out_file = open(path, 'w')

        for key, value in self.parameters.items():
            out_file.write('* {}: {}\n'.format(key, value))

        out_file.write('* data type: {}\n'.format(
            'mc' if self.is_sim else 'data'
        ))

        out_file.write('\n')

        for path in self.files:
            out_file.write(path)
            out_file.write('\n')

        out_file.close()

