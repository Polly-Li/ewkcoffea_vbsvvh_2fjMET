import setuptools

setuptools.setup(
    name='ewkcoffea',
    version='0.0.0',
    description='Analysis code for EWK analyses with coffea',
    packages=setuptools.find_packages(),
    # Include data files (Note: "include_package_data=True" does not seem to work)
    package_data={
        "ewkcoffea" : [
            "params/*",
        ],
    }
)

