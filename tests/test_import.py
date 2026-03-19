"""Simple import test for the package."""


def test_package_works():
    """Test that the package can be imported and basic functionality works."""
    # Test main package import
    import nemo_spinup_evaluation

    # Test version exists
    assert hasattr(nemo_spinup_evaluation, "__version__")
    assert nemo_spinup_evaluation.__version__ == "0.3.0"

    # Test CLI module can be imported
    from nemo_spinup_evaluation import cli

    assert hasattr(cli, "main")
    assert callable(cli.main)

    print("All imports successful!")
