## Pull Request

### Description
<!-- Provide a brief description of the changes in this PR -->

### Type of Change
<!-- Mark the relevant option with an 'x' -->

- [ ] Bug fix (non-breaking change which fixes an issue)
- [ ] New feature (non-breaking change which adds functionality)
- [ ] Breaking change (fix or feature that would cause existing functionality to not work as expected)
- [ ] Documentation update
- [ ] Refactor (code change that neither fixes a bug nor adds a feature)
- [ ] Performance improvement
- [ ] Configuration change
- [ ] Test addition or improvement

### Changes Made
<!-- List the specific changes made in this PR -->

-
-
-

### Testing
<!-- Describe the tests you ran to verify your changes -->

- [ ] Test suite passes (`pytest`)
- [ ] Application starts successfully (`python app.py`)
- [ ] Manual testing completed
- [ ] Docker build succeeds
- [ ] Single-analysis and batch flows both exercised (if the change touches either)

### Screenshots (if applicable)
<!-- Add screenshots to help explain your changes. Required for UI changes. -->

### Dependencies
<!-- List any new dependencies or changes to existing dependencies -->

- [ ] No new dependencies added
- [ ] Dependencies added (listed in `requirements.in` and compiled into `requirements.txt`)
- [ ] Dependencies updated

### Cross-repo impact
<!-- The analyser consumes the molAOP Builder's public REST API -->

- [ ] No impact on the Builder integration
- [ ] Changes the expected Builder API shape — `docs/KE-MAPPING-API-REFERENCE.md` updated
- [ ] Requires a matching change in `marvinm2/molAOP-builder` (link it below)

### Security Considerations

- [ ] No security implications
- [ ] Security review required
- [ ] New environment variables added
- [ ] Database changes made

### Checklist

- [ ] My code follows the project's coding standards
- [ ] I have performed a self-review of my code
- [ ] I have commented my code, particularly in hard-to-understand areas
- [ ] I have made corresponding changes to the documentation
- [ ] My changes generate no new warnings
- [ ] I have added tests that prove my fix is effective or that my feature works
- [ ] New and existing unit tests pass locally with my changes
- [ ] `CHANGELOG.md` updated for user-visible changes

### Related Issues
<!-- Link to any related issues.

     Use the KEYWORD, not a bare "(#123)". Only `Fixes`/`Closes`/`Resolves` auto-close
     the issue on merge; a bare reference just links to it. #122 stayed open through an
     entire release that shipped it because every commit referencing it used the bare
     form, and this repo has a longer history of the same: a past triage found five of
     thirteen open issues already delivered. If a PR completes an issue, say so with a
     keyword here. -->

Fixes #(issue number)
Closes #(issue number)
Related to #(issue number)

### Additional Notes
<!-- Add any additional notes for reviewers -->

---

**For Reviewers:**
- [ ] Code quality is acceptable
- [ ] Tests are comprehensive
- [ ] Documentation is updated
- [ ] Security considerations are addressed
- [ ] No obvious performance issues
