"""Simple import test for the package."""


def test_package_works():
    """Test that the package can be imported and basic functionality works."""
    # Test main package import
    import spinup_evaluation

    # Test version exists
    assert hasattr(spinup_evaluation, "__version__")
    assert spinup_evaluation.__version__ == "0.1.0"

    # Test CLI module can be imported
    from spinup_evaluation import cli

    assert hasattr(cli, "main")
    assert callable(cli.main)

    print("All imports successful!")
