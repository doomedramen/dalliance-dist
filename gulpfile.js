var gulp = require('gulp');
var rename = require('gulp-rename');
var uglify = require('gulp-uglify');
var fs = require('fs');

gulp.task('duplicate-js', function () {
  return gulp.src('./node_modules/dalliance/build/*.js')
    .pipe(gulp.dest('./dist'));
});

gulp.task('duplicate-css', function () {
  return gulp.src('./node_modules/dalliance/css/*')
    .pipe(gulp.dest('./dist/css'));
});

gulp.task('duplicate-fonts', function () {
  return gulp.src('./node_modules/dalliance/fonts/*')
    .pipe(gulp.dest('./dist/fonts'));
});

gulp.task('duplicate-img', function () {
  return gulp.src('./node_modules/dalliance/img/*')
    .pipe(gulp.dest('./dist/img'));
});

gulp.task('compress', ['duplicate'], function () {
  return gulp.src(['./dist/*.js', '!./dist/*.min.js'])
    .pipe(uglify())
    .pipe(rename({
      extname: '.min.js'
    }))
    .pipe(gulp.dest('dist'));
});

gulp.task('duplicate', ['duplicate-js', 'duplicate-css', 'duplicate-fonts', 'duplicate-img']);

gulp.task('update-bower', function () {
  var packageLocation = './node_modules/dalliance/package.json';
  var bowerLocation = './bower.json';
  var version = require(packageLocation).version;

  var file_content = fs.readFileSync(bowerLocation);
  var content = JSON.parse(file_content);
  content.version = version;
  fs.writeFileSync(bowerLocation, JSON.stringify(content, null, 2));
});

gulp.task('default', ['duplicate', 'compress', 'update-bower']);
