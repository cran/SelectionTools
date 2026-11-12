\name{rng.init}
\alias{rng.init}
\title{Initializes the random-number generators}
\description{Initializes the internal random-number generators, optionally
using a user-defined seed for reproducible runs.}
\usage{rng.init(seed = NULL)}
\arguments{
\item{seed}{\code{NULL} for ordinary initialization, or one finite integer
between \code{-.Machine$integer.max} and \code{.Machine$integer.max}.}
}
\details{
Package loading uses ordinary initialization, whose initial internal state is
clock-based. Without a supplied seed, the existing initialization sequence is
retained. Calling \code{rng.init()} after seeded use returns the internal
generators to clock-based initialization and restores the preceding OpenMP
dynamic-team setting; it does not reset R's random-number state.

With a seed, this function calls \code{set.seed(seed)}, draws 257 initialization
values from R's generator, and resets every internal generator, including its
buffered state. R's random-number state advances by these initialization draws.
One value seeds the original MathLib12 stream derivation, and the other values
seed the 256 internal simulation generators. R's currently selected RNG kind
is respected. Existing generators and simulation algorithms are unchanged.

Seeded mode disables OpenMP dynamic-team adjustment and uses static scheduling
for crossing, cross-validation and permutation tests. Reproducibility requires
the same inputs, sequence of calls, R RNG kind, build and numerical libraries,
and actual thread count. Thread configuration precedes initialization, and
OpenMP settings remain constant during the run. Normal top-level execution with up to
256 threads is supported; shared fallback execution is outside this guarantee.
A repeated complete experiment requires a separate reset of simulation data;
\code{reset.all()} does not reset random generators. Other functions with
explicit seed arguments may reset R's random-number state themselves.
}
\value{An invisible list returned by the compiled routine.
The function is primarily called for its initialization side effects.}
\examples{
st.set.num.threads(1)
rng.init(seed = 1001)
# A simulation runs after preparing its input data.
rng.init()
}
