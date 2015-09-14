var gulp = require('gulp');
var rename = require('gulp-rename');
var uglify = require('gulp-uglify');

gulp.task('duplicate', function () {
  return gulp.src('./node_modules/dalliance/build/*.js')
    .pipe(gulp.dest('./dist'));
});

gulp.task('compress', ['duplicate'], function () {
  return gulp.src(['./dist/*.js', '!./dist/*.min.js'])
    .pipe(uglify())
    .pipe(rename({
      extname: '.min.js'
    }))
    .pipe(gulp.dest('dist'));
});

gulp.task('default', ['duplicate', 'compress']);