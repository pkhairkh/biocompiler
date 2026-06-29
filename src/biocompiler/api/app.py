"""
BioCompiler API — FastAPI app factory.

Creates and configures the FastAPI application with middleware,
CORS, and route registration.
"""

import logging
import os
from typing import Any

from fastapi import FastAPI, Request
from fastapi.middleware.cors import CORSMiddleware

from . import auth as _auth_mod
from .models import MAX_REQUEST_SIZE
from .routes import _main_router, _protein_router

logger = logging.getLogger(__name__)


def get_cors_config() -> tuple[list[str], bool]:
    """Return (cors_origins, allow_credentials) from environment.

    This is used by the health endpoint and by create_app().
    """
    _raw_cors = os.environ.get("BIOCOMPILER_CORS_ORIGINS", "").strip()
    cors_origins = [o.strip() for o in _raw_cors.split(",") if o.strip()] if _raw_cors else []

    _raw_creds = os.environ.get("BIOCOMPILER_CORS_ALLOW_CREDENTIALS", "").strip().lower()
    _allow_creds = _raw_creds in ("true", "1", "yes") and cors_origins != ["*"]

    return cors_origins, _allow_creds


def create_app() -> FastAPI:
    """
    Create and configure the FastAPI application.

    Routes are registered from:
    - :mod:`biocompiler.api.routes._main_router` — core endpoints
    - :mod:`biocompiler.api.routes._protein_router` — protein analysis

    Middleware:
    - CORS (only when explicit origins are configured)
    - Request body size limit
    - Auth-mode warning header (optional mode)
    - Rate limiting
    """
    from .. import __version__

    cors_origins, _allow_creds = get_cors_config()

    app = FastAPI(
        title="BioCompiler API",
        description="Machine-verified gene design REST API with authentication and rate limiting",
        version=__version__,
        docs_url="/docs",
        redoc_url="/redoc",
    )

    # CORS middleware: only added when explicit origins are configured.
    if cors_origins:
        app.add_middleware(
            CORSMiddleware,
            allow_origins=cors_origins,
            allow_credentials=_allow_creds,
            allow_methods=["*"],
            allow_headers=["*"],
        )

    # Request body size middleware
    @app.middleware("http")
    async def request_size_limit_middleware(request: Request, call_next) -> Any:
        """Reject requests whose body exceeds MAX_REQUEST_SIZE bytes."""
        if request.method in ("POST", "PUT", "PATCH"):
            content_length = request.headers.get("content-length")
            if content_length is not None:
                try:
                    if int(content_length) > MAX_REQUEST_SIZE:
                        from fastapi.responses import JSONResponse
                        return JSONResponse(
                            status_code=413,
                            content={
                                "detail": (
                                    f"Request body too large ({int(content_length)} bytes). "
                                    f"Maximum: {MAX_REQUEST_SIZE} bytes."
                                ),
                            },
                        )
                except (ValueError, TypeError):
                    pass
        response = await call_next(request)
        return response

    # Auth-mode middleware — reads _AUTH_MODE from the auth module at runtime
    # so that test patches to _auth_mod._AUTH_MODE take effect immediately.
    @app.middleware("http")
    async def auth_mode_middleware(request: Request, call_next) -> Any:
        response = await call_next(request)
        if _auth_mod._AUTH_MODE == "optional":
            api_key = request.headers.get(_auth_mod.API_KEY_NAME)
            if not api_key or api_key not in _auth_mod._CONFIGURED_API_KEYS:
                response.headers["X-Auth-Warning"] = (
                    "Authentication is optional but recommended. "
                    "Unauthenticated access may be restricted in future versions."
                )
        return response

    # Rate limiting middleware
    @app.middleware("http")
    async def rate_limit_middleware(request: Request, call_next) -> Any:
        client_id = request.client.host if request.client else "unknown"
        try:
            _auth_mod._check_rate_limit(client_id)
        except Exception as e:
            from fastapi.responses import JSONResponse
            status_code = getattr(e, 'status_code', 429)
            detail = getattr(e, 'detail', str(e))
            return JSONResponse(status_code=status_code, content={"detail": detail})
        response = await call_next(request)
        return response

    # Register routers
    app.include_router(_main_router)
    app.include_router(_protein_router, prefix="/protein")

    return app
