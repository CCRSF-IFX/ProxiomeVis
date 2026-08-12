Review this entire repository as a senior R/Shiny engineer.

Your objectives:

1. Understand the application architecture
   - Identify entry points (app.R, ui.R/server.R, modules, golem structure, etc.)
   - Map data flow, reactive flow, and module dependencies
   - Summarize the app structure

2. Perform a comprehensive code audit
   - Bugs and potential runtime errors
   - Reactive programming issues
   - Performance bottlenecks
   - Memory inefficiencies
   - Security concerns
   - Maintainability issues
   - Scalability limitations
   - UX/UI problems

3. Evaluate Shiny best practices
   - Reactive design
   - Module organization
   - Separation of concerns
   - Naming conventions
   - Error handling
   - Logging
   - Testing coverage

4. Identify technical debt
   - Duplicated code
   - Overly complex functions
   - Dead code
   - Unused dependencies
   - Hard-coded values
   - Fragile assumptions

5. For each issue:
   - Explain why it is a problem
   - Estimate severity (Critical / High / Medium / Low)
   - Suggest a specific fix
   - Provide code examples when appropriate

6. Generate a prioritized improvement roadmap
   - Quick wins (<1 day)
   - Medium improvements (1–3 days)
   - Major refactors (>3 days)

7. Focus particularly on:
   - Shiny responsiveness
   - Large dataset handling
   - Async opportunities (future/promises)
   - Modularization
   - User experience
   - Deployment readiness

Before proposing changes:
- Read the entire repository.
- Build a mental model of the application.
- Do not modify code yet.
- First provide a detailed audit report.

Output:
A markdown report with:
1. Architecture overview
2. Major findings
3. Prioritized recommendations
4. Refactoring roadmap
5. Estimated impact of each recommendation
