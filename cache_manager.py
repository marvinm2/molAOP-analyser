"""
Cache manager for reference datasets and frequently accessed data.
"""
import time
import logging
from typing import Dict, Any, Optional, Callable
from threading import Lock

logger = logging.getLogger(__name__)

class CacheManager:
    """Simple in-memory cache with TTL support."""
    
    def __init__(self):
        self._cache: Dict[str, Dict[str, Any]] = {}
        self._lock = Lock()
    
    def get(self, key: str) -> Optional[Any]:
        """
        Get cached value by key.
        
        Args:
            key: Cache key
            
        Returns:
            Cached value or None if not found/expired
        """
        with self._lock:
            if key not in self._cache:
                return None
            
            cache_entry = self._cache[key]
            
            # Check if expired
            if cache_entry.get('ttl') and time.time() > cache_entry['expires_at']:
                del self._cache[key]
                logger.debug(f"Cache expired for key: {key}")
                return None
            
            logger.debug(f"Cache hit for key: {key}")
            return cache_entry['value']
    
    def set(self, key: str, value: Any, ttl: Optional[int] = None) -> None:
        """
        Set cached value.
        
        Args:
            key: Cache key
            value: Value to cache
            ttl: Time to live in seconds (None = no expiration)
        """
        with self._lock:
            cache_entry = {
                'value': value,
                'ttl': ttl,
                'created_at': time.time()
            }
            
            if ttl:
                cache_entry['expires_at'] = time.time() + ttl
            
            self._cache[key] = cache_entry
            logger.debug(f"Cache set for key: {key}, TTL: {ttl}")
    
    def get_or_set(self, key: str, factory: Callable[[], Any], ttl: Optional[int] = None) -> Any:
        """
        Get cached value or compute and cache it.
        
        Args:
            key: Cache key
            factory: Function to compute value if not cached
            ttl: Time to live in seconds
            
        Returns:
            Cached or computed value
        """
        value = self.get(key)
        
        if value is None:
            logger.debug(f"Cache miss for key: {key}, computing value")
            value = factory()
            self.set(key, value, ttl)
        
        return value
    
    def clear(self, pattern: Optional[str] = None) -> int:
        """
        Clear cache entries.
        
        Args:
            pattern: If provided, only clear keys containing this pattern
            
        Returns:
            Number of entries cleared
        """
        with self._lock:
            if pattern is None:
                count = len(self._cache)
                self._cache.clear()
                logger.info(f"Cleared entire cache: {count} entries")
                return count
            
            keys_to_remove = [k for k in self._cache.keys() if pattern in k]
            for key in keys_to_remove:
                del self._cache[key]
            
            logger.info(f"Cleared {len(keys_to_remove)} cache entries matching pattern: {pattern}")
            return len(keys_to_remove)
    
    def stats(self) -> Dict[str, Any]:
        """
        Get cache statistics.
        
        Returns:
            Dictionary with cache stats
        """
        with self._lock:
            total_entries = len(self._cache)
            expired_entries = 0
            current_time = time.time()
            
            for cache_entry in self._cache.values():
                if cache_entry.get('ttl') and current_time > cache_entry.get('expires_at', 0):
                    expired_entries += 1
            
            return {
                'total_entries': total_entries,
                'expired_entries': expired_entries,
                'active_entries': total_entries - expired_entries
            }

# Global cache instance
cache = CacheManager()

def cached_data_loader(cache_key: str, loader_func: Callable[[], Any], ttl: int = 3600) -> Any:
    """
    Decorator/utility for caching data loading operations.
    
    Args:
        cache_key: Unique key for this data
        loader_func: Function to load data if not cached
        ttl: Cache TTL in seconds (default 1 hour)
        
    Returns:
        Cached or loaded data
    """
    return cache.get_or_set(cache_key, loader_func, ttl)