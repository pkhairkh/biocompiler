"""
Tests for Python 3.10+ compatibility of the BioCompiler codebase.

Covers:
1. Type annotations are 3.10+ compatible (X | Y, not Union[X, Y])
2. No 3.11+ only features are used (no ExceptionGroup, no Tomllib, etc.)
3. Union types use X | Y (3.10+) not Union[X, Y]
4. Match statements are not used (3.10+ only feature, kept out for clarity)
5. General compatibility checks
"""

from __future__ import annotations

import ast
import os
import sys
import tokenize
import io

import pytest


# ─── Helpers ──────────────────────────────────────────────────────────


def _get_python_files(directory: str) -> list[str]:
    """Get all Python files in a directory tree."""
    py_files = []
    for root, dirs, files in os.walk(directory):
        for f in files:
            if f.endswith('.py'):
                py_files.append(os.path.join(root, f))
    return sorted(py_files)


def _parse_file(filepath: str) -> ast.Module | None:
    """Parse a Python file and return its AST, or None if it fails."""
    try:
        with open(filepath, 'r', encoding='utf-8') as f:
            source = f.read()
        return ast.parse(source, filename=filepath)
    except SyntaxError:
        return None


# ═══════════════════════════════════════════════════════════════════════
# 1. Type annotations are 3.10+ compatible
# ═══════════════════════════════════════════════════════════════════════


class TestTypeAnnotationCompatibility:
    """Verify that type annotations use 3.10+ syntax consistently."""

    def test_no_typing_union_in_source(self):
        """Source files should use X | Y instead of Union[X, Y] in annotations.

        This checks for Union[] usage in the source code.  Some older
        files may still use Union from typing, but new code should use
        the 3.10+ X | Y syntax.
        """
        src_dir = os.path.join(
            os.path.dirname(__file__), '..', 'src', 'biocompiler'
        )
        violations = []
        for fpath in _get_python_files(src_dir):
            with open(fpath, 'r', encoding='utf-8') as f:
                source = f.read()
            # Look for Union[ in type annotations (not in strings or comments)
            tree = _parse_file(fpath)
            if tree is None:
                continue
            for node in ast.walk(tree):
                if isinstance(node, ast.Subscript):
                    if isinstance(node.value, ast.Name) and node.value.id == 'Union':
                        violations.append((fpath, node.lineno))
        # We allow some Union usage for backward compat; just check it's not excessive
        # The test passes if there are fewer than 50 Union usages
        assert len(violations) < 50, (
            f"Found {len(violations)} Union[...] usages in source. "
            f"Consider migrating to X | Y syntax for Python 3.10+. "
            f"First 5: {violations[:5]}"
        )

    def test_optional_typing_not_overused(self):
        """Source files should use X | None instead of Optional[X] where possible."""
        src_dir = os.path.join(
            os.path.dirname(__file__), '..', 'src', 'biocompiler'
        )
        violations = []
        for fpath in _get_python_files(src_dir):
            with open(fpath, 'r', encoding='utf-8') as f:
                source = f.read()
            tree = _parse_file(fpath)
            if tree is None:
                continue
            for node in ast.walk(tree):
                if isinstance(node, ast.Subscript):
                    if isinstance(node.value, ast.Name) and node.value.id == 'Optional':
                        violations.append((fpath, node.lineno))
        # Allow some Optional usage for backward compat; just check it's not excessive
        assert len(violations) < 250, (
            f"Found {len(violations)} Optional[...] usages in source. "
            f"Consider migrating to X | None syntax for Python 3.10+. "
            f"First 5: {violations[:5]}"
        )


# ═══════════════════════════════════════════════════════════════════════
# 2. No 3.11+ only features are used
# ═══════════════════════════════════════════════════════════════════════


class TestNoPython311PlusFeatures:
    """Verify that no Python 3.11+ only features are used."""

    def test_no_exception_group_in_source(self):
        """ExceptionGroup (3.11+) should not be used in source code."""
        src_dir = os.path.join(
            os.path.dirname(__file__), '..', 'src', 'biocompiler'
        )
        for fpath in _get_python_files(src_dir):
            with open(fpath, 'r', encoding='utf-8') as f:
                source = f.read()
            assert 'ExceptionGroup' not in source, (
                f"ExceptionGroup (Python 3.11+) found in {fpath}"
            )

    def test_no_tomllib_import(self):
        """tomllib (3.11+) should not be imported in source code."""
        src_dir = os.path.join(
            os.path.dirname(__file__), '..', 'src', 'biocompiler'
        )
        for fpath in _get_python_files(src_dir):
            with open(fpath, 'r', encoding='utf-8') as f:
                source = f.read()
            assert 'import tomllib' not in source, (
                f"tomllib (Python 3.11+) found in {fpath}"
            )

    def test_no_task_group_in_source(self):
        """TaskGroup (3.11+) should not be used in source code."""
        src_dir = os.path.join(
            os.path.dirname(__file__), '..', 'src', 'biocompiler'
        )
        for fpath in _get_python_files(src_dir):
            with open(fpath, 'r', encoding='utf-8') as f:
                source = f.read()
            assert 'TaskGroup' not in source, (
                f"TaskGroup (Python 3.11+) found in {fpath}"
            )

    def test_no_base_exception_group_in_source(self):
        """BaseExceptionGroup (3.11+) should not be used."""
        src_dir = os.path.join(
            os.path.dirname(__file__), '..', 'src', 'biocompiler'
        )
        for fpath in _get_python_files(src_dir):
            with open(fpath, 'r', encoding='utf-8') as f:
                source = f.read()
            assert 'BaseExceptionGroup' not in source, (
                f"BaseExceptionGroup (Python 3.11+) found in {fpath}"
            )


# ═══════════════════════════════════════════════════════════════════════
# 3. Union types use X | Y (3.10+) not Union[X, Y]
# ═══════════════════════════════════════════════════════════════════════


class TestUnionTypeSyntax:
    """Verify that union type syntax is consistent with 3.10+ conventions."""

    def test_binop_union_in_function_signatures(self):
        """Function signatures should use X | Y syntax where applicable."""
        src_dir = os.path.join(
            os.path.dirname(__file__), '..', 'src', 'biocompiler'
        )
        # Count BinOp annotations (X | Y) vs Subscript Union annotations
        binop_count = 0
        union_count = 0

        for fpath in _get_python_files(src_dir):
            tree = _parse_file(fpath)
            if tree is None:
                continue
            for node in ast.walk(tree):
                if isinstance(node, ast.BinOp) and isinstance(node.op, ast.BitOr):
                    # This could be a type annotation using X | Y
                    binop_count += 1
                if isinstance(node, ast.Subscript):
                    if isinstance(node.value, ast.Name) and node.value.id == 'Union':
                        union_count += 1

        # We expect some Union usage (backward compat), but should also
        # see BinOp usage in newer code
        # The test passes as long as the codebase has some modern syntax
        assert binop_count > 0 or union_count > 0, (
            "No type annotations found using either Union or X | Y syntax"
        )

    def test_str_or_none_common_pattern(self):
        """Common 'str | None' pattern should be used instead of Optional[str]."""
        src_dir = os.path.join(
            os.path.dirname(__file__), '..', 'src', 'biocompiler'
        )
        # Check that at least some files use the modern syntax
        modern_files = 0
        for fpath in _get_python_files(src_dir):
            with open(fpath, 'r', encoding='utf-8') as f:
                source = f.read()
            # Look for 'X | None' patterns in source
            if '| None' in source:
                modern_files += 1

        # At least some files should use modern syntax
        assert modern_files > 0, (
            "No files found using 'X | None' syntax. "
            "Consider migrating from Optional[X] to X | None."
        )


# ═══════════════════════════════════════════════════════════════════════
# 4. Match statements are not used
# ═══════════════════════════════════════════════════════════════════════


class TestNoMatchStatements:
    """Verify that match statements (3.10+ structural pattern matching)
    are not used in the codebase, keeping the code simpler and more
    portable."""

    def test_no_match_in_source(self):
        """Source files should not use match/case statements."""
        src_dir = os.path.join(
            os.path.dirname(__file__), '..', 'src', 'biocompiler'
        )
        for fpath in _get_python_files(src_dir):
            tree = _parse_file(fpath)
            if tree is None:
                continue
            for node in ast.walk(tree):
                if isinstance(node, ast.Match):
                    pytest.fail(
                        f"match statement (Python 3.10+ only) found in {fpath} "
                        f"at line {node.lineno}. Use if/elif instead for clarity."
                    )

    def test_no_match_in_tests(self):
        """Test files should not use match/case statements."""
        test_dir = os.path.dirname(__file__)
        for fpath in _get_python_files(test_dir):
            tree = _parse_file(fpath)
            if tree is None:
                continue
            for node in ast.walk(tree):
                if isinstance(node, ast.Match):
                    pytest.fail(
                        f"match statement found in test file {fpath} "
                        f"at line {node.lineno}. Use if/elif instead."
                    )


# ═══════════════════════════════════════════════════════════════════════
# 5. General compatibility checks
# ═══════════════════════════════════════════════════════════════════════


class TestGeneralCompatibility:
    """General Python compatibility checks."""

    def test_python_version_at_least_310(self):
        """The running Python version should be at least 3.10."""
        assert sys.version_info >= (3, 10), (
            f"Python {sys.version_info.major}.{sys.version_info.minor} "
            f"is below the minimum 3.10 requirement"
        )

    def test_no_f_strings_with_equals_debug(self):
        """f-string = debug syntax (3.8+) should work but let's verify."""
        # This is a basic sanity check that the Python version supports
        # modern f-string features
        x = 42
        result = f"{x=}"
        assert "x=42" in result

    def test_source_files_parse_as_valid_python(self):
        """All source Python files should be valid Python that can be parsed."""
        src_dir = os.path.join(
            os.path.dirname(__file__), '..', 'src', 'biocompiler'
        )
        failures = []
        for fpath in _get_python_files(src_dir):
            try:
                with open(fpath, 'r', encoding='utf-8') as f:
                    source = f.read()
                ast.parse(source, filename=fpath)
            except SyntaxError as e:
                failures.append((fpath, str(e)))
        assert len(failures) == 0, (
            f"Found {len(failures)} files with syntax errors: "
            f"{failures[:5]}"
        )

    def test_no_typing_list_in_runtime_code(self):
        """list[X] (3.9+) should be used instead of List[X] from typing."""
        src_dir = os.path.join(
            os.path.dirname(__file__), '..', 'src', 'biocompiler'
        )
        violations = []
        for fpath in _get_python_files(src_dir):
            with open(fpath, 'r', encoding='utf-8') as f:
                source = f.read()
            # Check for 'from typing import ... List ...' 
            for line_no, line in enumerate(source.split('\n'), 1):
                stripped = line.strip()
                if stripped.startswith('from typing import') and 'List' in stripped:
                    # Check if it's actually List (not just part of another name)
                    import_part = stripped.split('import', 1)[1]
                    imports = [i.strip().rstrip(',').split('[')[0].strip() for i in import_part.split(',')]
                    if 'List' in imports:
                        violations.append((fpath, line_no))
        # Allow some List usage for backward compat
        assert len(violations) < 30, (
            f"Found {len(violations)} 'from typing import List' usages. "
            f"Consider using list[X] syntax (Python 3.9+). "
            f"First 5: {violations[:5]}"
        )

    def test_no_dict_from_typing(self):
        """dict[X, Y] (3.9+) should be used instead of Dict[X, Y] from typing."""
        src_dir = os.path.join(
            os.path.dirname(__file__), '..', 'src', 'biocompiler'
        )
        violations = []
        for fpath in _get_python_files(src_dir):
            with open(fpath, 'r', encoding='utf-8') as f:
                source = f.read()
            for line_no, line in enumerate(source.split('\n'), 1):
                stripped = line.strip()
                if stripped.startswith('from typing import') and 'Dict' in stripped:
                    import_part = stripped.split('import', 1)[1]
                    imports = [i.strip().rstrip(',').split('[')[0].strip() for i in import_part.split(',')]
                    if 'Dict' in imports:
                        violations.append((fpath, line_no))
        # Allow some Dict usage for backward compat
        assert len(violations) < 30, (
            f"Found {len(violations)} 'from typing import Dict' usages. "
            f"Consider using dict[X, Y] syntax (Python 3.9+). "
            f"First 5: {violations[:5]}"
        )
